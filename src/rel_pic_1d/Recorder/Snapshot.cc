/*
 * Copyright (c) 2020-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "Snapshot.h"

#include <algorithm>
#include <functional>
#include <iterator>
#include <stdexcept>
#include <type_traits>

namespace {
template <class Tuple> struct Hash;
template <class T> Hash(T const &t) -> Hash<T>;
//
template <class... Ts> struct Hash<std::tuple<Ts...>> {
    std::tuple<Ts...> const t;

    [[nodiscard]] constexpr std::size_t operator()() const noexcept
    {
        return hash(std::index_sequence_for<Ts...>{});
    }

private:
    template <std::size_t... Is>
    [[nodiscard]] constexpr std::size_t hash(std::index_sequence<Is...>) const noexcept
    {
        std::size_t hash = 0;
        return (..., ((hash <<= 1) ^= this->hash(std::get<Is>(t))));
    }
    template <class T> [[nodiscard]] static constexpr std::size_t hash(T const &x) noexcept
    {
        return std::hash<T>{}(x);
    }
};
} // namespace

// MARK:- P1D::Snapshot
//
P1D::Snapshot::Snapshot(parallel::mpi::Comm _comm, ParamSet const &params)
: comm{ std::move(_comm), tag }
, signature{ Hash{ serialize(params) }() }
, wd{ params.working_directory }
{
    if (!comm->operator bool())
        throw std::invalid_argument{ std::string{ __PRETTY_FUNCTION__ } + " - invalid mpi::Comm" };

    // method dispatch
    //
    if (is_master()) {
        save = &Snapshot::save_master;
        load = &Snapshot::load_master;
    } else {
        save = &Snapshot::save_worker;
        load = &Snapshot::load_worker;
    }
}

inline std::string P1D::Snapshot::filepath() const
{
    constexpr char filename[] = "snapshot.h5";
    return wd + "/" + filename;
}

template <class MType, long N>
auto P1D::Snapshot::save_helper(hdf5::Group &root, GridQ<MType, N> const &grid,
                                std::string const &basename) const -> hdf5::Dataset
{
    using FType = Real;
    static_assert(alignof(MType) == alignof(FType), "memory and file type mis-alignment");
    static_assert(0 == sizeof(MType) % sizeof(FType), "memory and file type size incompatible");
    constexpr auto len   = sizeof(MType) / sizeof(FType);
    auto const     ftype = hdf5::make_type<FType>();

    auto payload = *comm.gather<MType>({ grid.begin(), grid.end() }, master);
    auto mspace  = hdf5::Space::simple({ payload.size(), len });

    // fixed dataset
    auto fspace = hdf5::Space::simple({ payload.size(), len });
    auto dset   = root.dataset(basename.c_str(), ftype, fspace);

    // export
    mspace.select_all();
    fspace.select_all();
    dset.write(fspace, payload.data(), ftype, mspace);

    return dset;
}
void P1D::Snapshot::save_helper(hdf5::Group &root, PartSpecies const &sp) const
{
    // collect
    auto payload = sp.dump_ptls();
    auto tk      = comm.ibsend(std::move(payload), master);
    payload.clear();
    payload.reserve(static_cast<unsigned long>(sp->Nc * sp.params.Nx));
    for (int rank = 0, size = comm.size(); rank < size; ++rank) {
        comm.recv<Particle>({}, rank).unpack(
            [](auto incoming, auto &payload) {
                payload.insert(payload.end(), std::make_move_iterator(begin(incoming)),
                               std::make_move_iterator(end(incoming)));
            },
            payload);
    }
    std::move(tk).wait();

    // export
    static_assert(sizeof(Particle) % sizeof(Real) == 0);
    static_assert(alignof(Particle) == alignof(Real));
    auto const real_type = hdf5::make_type<Real>();
    auto mspace = hdf5::Space::simple({ payload.size(), sizeof(Particle) / real_type.size() });
    {
        mspace.select(H5S_SELECT_SET, { 0U, 0U }, { payload.size(), 3U });
        auto fspace = hdf5::Space::simple({ payload.size(), 3U });
        auto dset = root.dataset("vel", real_type, fspace, hdf5::PList::dapl(), hdf5::PList::dcpl())
                 << sp;
        fspace.select_all();
        dset.write(fspace, payload.data(), real_type, mspace);
    }
    {
        mspace.select(H5S_SELECT_SET, { 0U, 3U }, { payload.size(), 1U });
        auto fspace = hdf5::Space::simple({ payload.size(), 1U });
        auto dset
            = root.dataset("pos_x", real_type, fspace, hdf5::PList::dapl(), hdf5::PList::dcpl())
           << sp;
        fspace.select_all();
        dset.write(fspace, payload.data(), real_type, mspace);
    }
    {
        mspace.select(H5S_SELECT_SET, { 0U, 4U }, { payload.size(), 2U });
        auto fspace = hdf5::Space::simple({ payload.size(), 2U });
        auto dset = root.dataset("fw", real_type, fspace, hdf5::PList::dapl(), hdf5::PList::dcpl())
                 << sp;
        fspace.select_all();
        dset.write(fspace, payload.data(), real_type, mspace);
    }
    {
        mspace.select(H5S_SELECT_SET, { 0U, 6U }, { payload.size(), 1U });
        auto fspace = hdf5::Space::simple(payload.size());
        auto dset
            = root.dataset("gamma", real_type, fspace, hdf5::PList::dapl(), hdf5::PList::dcpl())
           << sp;
        fspace.select_all();
        dset.write(fspace, payload.data(), real_type, mspace);
    }
}
void P1D::Snapshot::save_master(Domain const &domain, long const step_count) const &
{
    // create hdf5 file and root group
    hdf5::Group root = hdf5::File(hdf5::File::trunc_tag{}, filepath().c_str())
                           .group("pic_1d", hdf5::PList::gapl(), hdf5::PList::gcpl());
    root << domain.params;

    // step_count & signature
    root.attribute("step_count", hdf5::make_type(step_count), hdf5::Space::scalar())
        .write(step_count);
    root.attribute("signature", hdf5::make_type(signature), hdf5::Space::scalar()).write(signature);

    // B & E
    save_helper(root, domain.bfield, "bfield") << domain.bfield;
    save_helper(root, domain.efield, "efield") << domain.efield;

    // particles
    for (unsigned i = 0; i < domain.part_species.size(); ++i) {
        PartSpecies const &sp    = domain.part_species.at(i);
        std::string const  gname = std::string{ "part_species[" } + std::to_string(i) + "]";
        auto group = root.group(gname.c_str(), hdf5::PList::gapl(), hdf5::PList::gcpl()) << sp;
        save_helper(group, sp);
    }

    // cold fluid
    for (unsigned i = 0; i < domain.cold_species.size(); ++i) {
        ColdSpecies const &sp    = domain.cold_species.at(i);
        std::string const  gname = std::string{ "cold_species[" } + std::to_string(i) + "]";
        auto group = root.group(gname.c_str(), hdf5::PList::gapl(), hdf5::PList::gcpl()) << sp;
        save_helper(group, sp.mom0_full, "mom0_full") << sp;
        save_helper(group, sp.mom1_full, "mom1_full") << sp;
    }

    root.flush();
}
void P1D::Snapshot::save_worker(Domain const &domain, long) const &
{
    // B & E
    comm.gather<1>(domain.bfield.begin(), domain.bfield.end(), nullptr, master);
    comm.gather<1>(domain.efield.begin(), domain.efield.end(), nullptr, master);

    // particles
    for (PartSpecies const &sp : domain.part_species) {
        comm.ibsend(sp.dump_ptls(), master).wait();
    }

    // cold fluid
    for (ColdSpecies const &sp : domain.cold_species) {
        comm.gather<0>(sp.mom0_full.begin(), sp.mom0_full.end(), nullptr, master);
        comm.gather<1>(sp.mom1_full.begin(), sp.mom1_full.end(), nullptr, master);
    }
}

template <class MType, long N>
auto P1D::Snapshot::load_helper(hdf5::Group const &root, GridQ<MType, N> &grid,
                                std::string const &basename) const -> hdf5::Dataset
{
    using FType = Real;
    static_assert(alignof(MType) == alignof(FType), "memory and file type mis-alignment");
    static_assert(0 == sizeof(MType) % sizeof(FType), "memory and file type size incompatible");
    constexpr auto len   = sizeof(MType) / sizeof(FType);
    auto const     ftype = hdf5::make_type<FType>();

    std::vector<MType> payload(static_cast<unsigned long>(grid.size() * comm.size()));
    auto               mspace = hdf5::Space::simple({ payload.size(), len });

    // open dataset
    auto       dset   = root.dataset(basename.c_str());
    auto       fspace = dset.space();
    auto const extent = fspace.simple_extent().first;
    if (extent.rank() != 2 || extent[0] != payload.size() || extent[1] != len)
        throw std::runtime_error{ std::string{ __PRETTY_FUNCTION__ }
                                  + " - incompatible extent : " + basename };

    // import
    mspace.select_all();
    fspace.select_all();
    dset.read(fspace, payload.data(), ftype, mspace);

    // distribute
    comm.scatter(payload.data(), grid.begin(), grid.end(), master);

    return dset;
}
void P1D::Snapshot::load_helper(hdf5::Group const &root, PartSpecies &sp) const
{
    std::vector<Particle> payload(static_cast<unsigned long>(sp->Nc * sp.params.Nx));

    // import
    static_assert(sizeof(Particle) % sizeof(Real) == 0);
    static_assert(alignof(Particle) == alignof(Real));
    auto const real_type = hdf5::make_type<Real>();
    auto mspace = hdf5::Space::simple({ payload.size(), sizeof(Particle) / real_type.size() });
    {
        auto       dset   = root.dataset("vel");
        auto       fspace = dset.space();
        auto const extent = fspace.simple_extent().first;
        if (extent.rank() != 2 || extent[0] != payload.size() || extent[1] != 3)
            throw std::runtime_error{ std::string{ __PRETTY_FUNCTION__ }
                                      + " - incompatible extent : vel" };

        mspace.select(H5S_SELECT_SET, { 0U, 0U }, { payload.size(), 3U });
        fspace.select_all();
        dset.read(fspace, payload.data(), real_type, mspace);
    }
    {
        auto       dset   = root.dataset("pos_x");
        auto       fspace = dset.space();
        auto const extent = fspace.simple_extent().first;
        if (extent.rank() != 2 || extent[0] != payload.size() || extent[1] != 1)
            throw std::runtime_error{ std::string{ __PRETTY_FUNCTION__ }
                                      + " - incompatible extent : pos_x" };

        mspace.select(H5S_SELECT_SET, { 0U, 3U }, { payload.size(), 1U });
        fspace.select_all();
        dset.read(fspace, payload.data(), real_type, mspace);
    }
    {
        auto       dset   = root.dataset("fw");
        auto       fspace = dset.space();
        auto const extent = fspace.simple_extent().first;
        if (extent.rank() != 2 || extent[0] != payload.size() || extent[1] != 2)
            throw std::runtime_error{ std::string{ __PRETTY_FUNCTION__ }
                                      + " - incompatible extent : fw" };

        mspace.select(H5S_SELECT_SET, { 0U, 4U }, { payload.size(), 2U });
        fspace.select_all();
        dset.read(fspace, payload.data(), real_type, mspace);
    }
    {
        auto       dset   = root.dataset("gamma");
        auto       fspace = dset.space();
        auto const extent = fspace.simple_extent().first;
        if (extent.rank() != 1 || extent[0] != payload.size())
            throw std::runtime_error{ std::string{ __PRETTY_FUNCTION__ }
                                      + " - incompatible extent : gamma" };

        mspace.select(H5S_SELECT_SET, { 0U, 6U }, { payload.size(), 1U });
        fspace.select_all();
        dset.read(fspace, payload.data(), real_type, mspace);
    }

    // distribute
    std::reverse(payload.begin(), payload.end()); // This is to make the sequence the same as in
                                                  // the version of multi-thread particle loading.
    auto last = payload.crbegin();
    for (long i = 0; i < sp.params.Nx; ++i) { // assumption here is that the number of particles is
                                              // divisible to the number of grid points.
        std::advance(last, sp->Nc);
        comm.bcast<Particle>({ payload.crbegin(), last }, master)
            .unpack(
                [](auto payload, PartSpecies &sp, bool append) {
                    sp.load_ptls(std::move(payload), append);
                },
                sp, i);
        payload.erase(last.base(), payload.end());
    }
    if (!payload.empty())
        throw std::runtime_error{ std::string{ __PRETTY_FUNCTION__ }
                                  + " - particles still remaining after distribution" };
}
long P1D::Snapshot::load_master(Domain &domain) const &
{
    // open hdf5 file and root group
    hdf5::Group root = hdf5::File(hdf5::File::rdonly_tag{}, filepath().c_str()).group("pic_1d");

    // verify signature
    if (signature != root.attribute("signature").read<decltype(signature)>())
        throw std::runtime_error{ std::string{ __PRETTY_FUNCTION__ }
                                  + " - signature verification failed" };

    // B & E
    load_helper(root, domain.bfield, "bfield");
    load_helper(root, domain.efield, "efield");

    // particles
    for (unsigned i = 0; i < domain.part_species.size(); ++i) {
        PartSpecies      &sp    = domain.part_species.at(i);
        std::string const gname = std::string{ "part_species[" } + std::to_string(i) + "]";
        auto const        group = root.group(gname.c_str());
        load_helper(group, sp);

        auto const count = *comm.all_reduce<long>(parallel::mpi::ReduceOp::plus(),
                                                  static_cast<long>(sp.bucket.size()));
        if (sp->Nc * sp.params.Nx != count)
            throw std::runtime_error{ std::string{ __PRETTY_FUNCTION__ }
                                      + " - particle count inconsistent for species "
                                      + std::to_string(i) };
    }

    // cold fluid
    for (unsigned i = 0; i < domain.cold_species.size(); ++i) {
        ColdSpecies      &sp    = domain.cold_species.at(i);
        std::string const gname = std::string{ "cold_species[" } + std::to_string(i) + "]";
        auto const        group = root.group(gname.c_str());
        load_helper(group, sp.mom0_full, "mom0_full");
        load_helper(group, sp.mom1_full, "mom1_full");
    }

    // step count
    auto const step_count = root.attribute("step_count").read<long>();
    return comm.bcast<long>(step_count, master).unpack([step_count](auto) {
        return step_count;
    });
}
long P1D::Snapshot::load_worker(Domain &domain) const &
{
    // B & E
    comm.scatter<1>(nullptr, domain.bfield.begin(), domain.bfield.end(), master);
    comm.scatter<1>(nullptr, domain.efield.begin(), domain.efield.end(), master);

    // particles
    for (PartSpecies &sp : domain.part_species) {
        std::vector<Particle> payload(static_cast<unsigned long>(sp->Nc));
        for (long i = 0; i < sp.params.Nx; ++i) {
            comm.bcast<Particle>(payload, master)
                .unpack(
                    [](auto payload, PartSpecies &sp, bool append) {
                        sp.load_ptls(std::move(payload), append);
                    },
                    sp, i);
        }
        comm.all_reduce<long>(parallel::mpi::ReduceOp::plus(), static_cast<long>(sp.bucket.size()))
            .unpack([](auto) {});
    }

    // cold fluid
    for (ColdSpecies &sp : domain.cold_species) {
        comm.scatter<0>(nullptr, sp.mom0_full.begin(), sp.mom0_full.end(), master);
        comm.scatter<1>(nullptr, sp.mom1_full.begin(), sp.mom1_full.end(), master);
    }

    // step count
    return comm.bcast<long>({}, master);
}
