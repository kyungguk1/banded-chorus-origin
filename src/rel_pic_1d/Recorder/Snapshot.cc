/*
 * Copyright (c) 2020-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "Snapshot.h"

#include <algorithm>
#include <cstddef> // offsetof
#include <functional>
#include <iterator>
#include <stdexcept>
#include <type_traits>

PIC1D_BEGIN_NAMESPACE
namespace {
template <class Tuple>
struct Hash;
template <class T>
Hash(T const &t) -> Hash<T>;
//
template <class... Ts>
struct Hash<std::tuple<Ts...>> {
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
    template <class T>
    [[nodiscard]] static constexpr std::size_t hash(T const &x) noexcept
    {
        return std::hash<T>{}(x);
    }
};
} // namespace

Snapshot::Snapshot(parallel::mpi::Comm _comm, ParamSet const &params, long const subdomain_color)
: comm{ std::move(_comm) }
, signature{ Hash{ serialize(params) }() }
, wd{ params.working_directory }
{
    if (!comm->operator bool())
        throw std::invalid_argument{ std::string{ __PRETTY_FUNCTION__ } + " - invalid mpi::Comm" };

    if (subdomain_color >= 0)
        (filename_suffix = '_') += std::to_string(subdomain_color);

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

inline std::string Snapshot::filepath() const
{
    constexpr char basename[]  = "snapshot";
    constexpr char extension[] = ".h5";
    return wd + "/" + basename + filename_suffix + extension;
}

template <class T, long N>
auto Snapshot::save_helper(hdf5::Group &root, Grid<T, N, Pad> const &grid, std::string const &basename) const -> hdf5::Dataset
{
    static_assert(alignof(T) == alignof(Real), "memory and file type mis-alignment");
    static_assert(0 == sizeof(T) % sizeof(Real), "memory and file type size incompatible");
    constexpr auto len  = sizeof(T) / sizeof(Real);
    auto const     type = hdf5::make_type<Real>();

    std::vector<T> payload = *comm.gather<T>({ grid.begin(), grid.end() }, master);
    auto           mspace  = hdf5::Space::simple({ payload.size(), len });

    // fixed dataset
    auto fspace = hdf5::Space::simple({ payload.size(), len });
    auto dset   = root.dataset(basename.c_str(), type, fspace);

    // export
    mspace.select_all();
    fspace.select_all();
    dset.write(fspace, payload.data(), type, mspace);

    return dset;
}
void Snapshot::save_helper(hdf5::Group &root, PartSpecies const &sp) const
{
    // collect
    std::vector<Particle> payload;
    long const            Np = sp->Nc * sp.params.Nx / sp.params.number_of_distributed_particle_subdomain_clones;
    payload.reserve(static_cast<unsigned long>(Np));
    auto tk = comm.ibsend(sp.dump_ptls(), { master, tag });
    for (int rank = 0, size = comm.size(); rank < size; ++rank) {
        comm.recv<Particle>({}, { rank, tag })
            .unpack(
                [](auto incoming, auto &payload) {
                    payload.insert(payload.end(), std::make_move_iterator(begin(incoming)),
                                   std::make_move_iterator(end(incoming)));
                },
                payload);
    }
    std::move(tk).wait();

    // export
    constexpr auto unit_size = sizeof(Real);
    static_assert(sizeof(Particle) % unit_size == 0);
    static_assert(alignof(Particle) == alignof(Real));
    auto mspace = hdf5::Space::simple({ payload.size(), sizeof(Particle) / unit_size });
    {
        auto const type = hdf5::make_type<Real>();
        using T         = std::decay_t<decltype(std::declval<Particle>().g_vel)>;

        hdf5::Extent const start = { 0U, offsetof(Particle, g_vel) / unit_size };
        hdf5::Extent const count = { payload.size(), sizeof(T) / unit_size };
        mspace.select(H5S_SELECT_SET, start, count);
        auto fspace = hdf5::Space::simple(count);

        auto dset = root.dataset("vel", type, fspace, hdf5::PList::dapl(), hdf5::PList::dcpl()) << sp;
        fspace.select_all();
        dset.write(fspace, payload.data(), type, mspace);
    }
    {
        auto const type = hdf5::make_type<Real>();
        using T         = std::decay_t<decltype(std::declval<Particle>().pos)>;

        hdf5::Extent const start = { 0U, offsetof(Particle, pos) / unit_size };
        hdf5::Extent const count = { payload.size(), sizeof(T) / unit_size };
        mspace.select(H5S_SELECT_SET, start, count);
        auto fspace = hdf5::Space::simple(count);

        auto dset = root.dataset("pos", type, fspace, hdf5::PList::dapl(), hdf5::PList::dcpl()) << sp;
        fspace.select_all();
        dset.write(fspace, payload.data(), type, mspace);
    }
    {
        auto const type = hdf5::make_type<Real>();
        using T         = std::decay_t<decltype(std::declval<Particle>().psd)>;

        hdf5::Extent const start = { 0U, offsetof(Particle, psd) / unit_size };
        hdf5::Extent const count = { payload.size(), sizeof(T) / unit_size };
        mspace.select(H5S_SELECT_SET, start, count);
        auto fspace = hdf5::Space::simple(count);

        auto dset = root.dataset("psd", type, fspace, hdf5::PList::dapl(), hdf5::PList::dcpl()) << sp;
        fspace.select_all();
        dset.write(fspace, payload.data(), type, mspace);
    }
    {
        auto const type = hdf5::make_type<Real>();
        using T         = std::decay_t<decltype(std::declval<Particle>().gamma)>;

        hdf5::Extent const start = { 0U, offsetof(Particle, gamma) / unit_size };
        hdf5::Extent const count = { payload.size(), sizeof(T) / unit_size };
        mspace.select(H5S_SELECT_SET, start, count);
        auto fspace = hdf5::Space::simple(payload.size());

        auto dset = root.dataset("gamma", type, fspace, hdf5::PList::dapl(), hdf5::PList::dcpl()) << sp;
        fspace.select_all();
        dset.write(fspace, payload.data(), type, mspace);
    }
    {
        auto const type = hdf5::make_type<long>();
        using T         = std::decay_t<decltype(std::declval<Particle>().id)>;

        hdf5::Extent const start = { 0U, offsetof(Particle, id) / unit_size };
        hdf5::Extent const count = { payload.size(), sizeof(T) / unit_size };
        mspace.select(H5S_SELECT_SET, start, count);
        auto fspace = hdf5::Space::simple(payload.size());

        auto dset = root.dataset("id", type, fspace, hdf5::PList::dapl(), hdf5::PList::dcpl()) << sp;
        fspace.select_all();
        dset.write(fspace, payload.data(), type, mspace);
    }
}
void Snapshot::save_master(Domain const &domain, long const step_count) const &
{
    // create hdf5 file and root group
    hdf5::Group root = hdf5::File(hdf5::File::trunc_tag{}, filepath().c_str())
                           .group("rel_pic_1d", hdf5::PList::gapl(), hdf5::PList::gcpl());
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
        std::string const  gname = std::string{ "part_species" } + '@' + std::to_string(i);
        auto               group = root.group(gname.c_str(), hdf5::PList::gapl(), hdf5::PList::gcpl()) << sp;
        save_helper(group, sp);
    }

    // cold fluid
    for (unsigned i = 0; i < domain.cold_species.size(); ++i) {
        ColdSpecies const &sp    = domain.cold_species.at(i);
        std::string const  gname = std::string{ "cold_species" } + '@' + std::to_string(i);
        auto               group = root.group(gname.c_str(), hdf5::PList::gapl(), hdf5::PList::gcpl()) << sp;
        save_helper(group, sp.mom0_full, "mom0_full") << sp;
        save_helper(group, sp.mom1_full, "mom1_full") << sp;
    }

    root.flush();
}
void Snapshot::save_worker(Domain const &domain, long) const &
{
    // B & E
    comm.gather<1>(domain.bfield.begin(), domain.bfield.end(), nullptr, master);
    comm.gather<1>(domain.efield.begin(), domain.efield.end(), nullptr, master);

    // particles
    for (PartSpecies const &sp : domain.part_species) {
        comm.ibsend(sp.dump_ptls(), { master, tag }).wait();
    }

    // cold fluid
    for (ColdSpecies const &sp : domain.cold_species) {
        comm.gather<0>(sp.mom0_full.begin(), sp.mom0_full.end(), nullptr, master);
        comm.gather<1>(sp.mom1_full.begin(), sp.mom1_full.end(), nullptr, master);
    }
}

template <class T, long N>
void Snapshot::load_helper(hdf5::Group const &root, Grid<T, N, Pad> &grid, std::string const &basename) const
{
    static_assert(alignof(T) == alignof(Real), "memory and file type mis-alignment");
    static_assert(0 == sizeof(T) % sizeof(Real), "memory and file type size incompatible");
    constexpr auto len  = sizeof(T) / sizeof(Real);
    auto const     type = hdf5::make_type<Real>();

    std::vector<T> payload(static_cast<unsigned long>(grid.size() * comm.size()));
    auto           mspace = hdf5::Space::simple({ payload.size(), len });

    // open dataset
    auto       dset   = root.dataset(basename.c_str());
    auto       fspace = dset.space();
    auto const extent = fspace.simple_extent().first;
    if (extent.rank() != 2 || extent[0] != payload.size() || extent[1] != len)
        throw std::runtime_error{ std::string{ __PRETTY_FUNCTION__ } + " - incompatible extent : " + basename };

    // import
    mspace.select_all();
    fspace.select_all();
    dset.read(fspace, payload.data(), type, mspace);

    // distribute
    comm.scatter(payload.data(), grid.begin(), grid.end(), master);
}
void Snapshot::load_helper(hdf5::Group const &root, PartSpecies &sp) const
{
    std::vector<Particle> payload;

    // get particle count
    {
        auto const extent = root.dataset("id").space().simple_extent().first;
        if (extent.rank() != 1)
            throw std::runtime_error{ std::string{ __PRETTY_FUNCTION__ } + " - expects 1D array in obtaining particle count" };
        payload.resize(extent[0]);
    }

    // import
    constexpr auto unit_size = sizeof(Real);
    static_assert(sizeof(Particle) % unit_size == 0);
    static_assert(alignof(Particle) == alignof(Real));
    auto mspace = hdf5::Space::simple({ payload.size(), sizeof(Particle) / unit_size });
    {
        auto const type   = hdf5::make_type<Real>();
        using T           = std::decay_t<decltype(std::declval<Particle>().g_vel)>;
        auto       dset   = root.dataset("vel");
        auto       fspace = dset.space();
        auto const extent = fspace.simple_extent().first;
        if (extent.rank() != 2 || extent[0] != payload.size() || extent[1] != sizeof(T) / unit_size)
            throw std::runtime_error{ std::string{ __PRETTY_FUNCTION__ } + " - incompatible extent : vel" };

        hdf5::Extent const start = { 0U, offsetof(Particle, g_vel) / unit_size };
        hdf5::Extent const count = { payload.size(), sizeof(T) / unit_size };
        mspace.select(H5S_SELECT_SET, start, count);
        fspace.select_all();
        dset.read(fspace, payload.data(), type, mspace);
    }
    {
        auto const type   = hdf5::make_type<Real>();
        using T           = std::decay_t<decltype(std::declval<Particle>().pos)>;
        auto       dset   = root.dataset("pos");
        auto       fspace = dset.space();
        auto const extent = fspace.simple_extent().first;
        if (extent.rank() != 2 || extent[0] != payload.size() || extent[1] != sizeof(T) / unit_size)
            throw std::runtime_error{ std::string{ __PRETTY_FUNCTION__ } + " - incompatible extent : pos" };

        hdf5::Extent const start = { 0U, offsetof(Particle, pos) / unit_size };
        hdf5::Extent const count = { payload.size(), sizeof(T) / unit_size };
        mspace.select(H5S_SELECT_SET, start, count);
        fspace.select_all();
        dset.read(fspace, payload.data(), type, mspace);
    }
    {
        auto const type   = hdf5::make_type<Real>();
        using T           = std::decay_t<decltype(std::declval<Particle>().psd)>;
        auto       dset   = root.dataset("psd");
        auto       fspace = dset.space();
        auto const extent = fspace.simple_extent().first;
        if (extent.rank() != 2 || extent[0] != payload.size() || extent[1] != sizeof(T) / unit_size)
            throw std::runtime_error{ std::string{ __PRETTY_FUNCTION__ } + " - incompatible extent : psd" };

        hdf5::Extent const start = { 0U, offsetof(Particle, psd) / unit_size };
        hdf5::Extent const count = { payload.size(), sizeof(T) / unit_size };
        mspace.select(H5S_SELECT_SET, start, count);
        fspace.select_all();
        dset.read(fspace, payload.data(), type, mspace);
    }
    {
        auto const type   = hdf5::make_type<Real>();
        using T           = std::decay_t<decltype(std::declval<Particle>().gamma)>;
        auto       dset   = root.dataset("gamma");
        auto       fspace = dset.space();
        auto const extent = fspace.simple_extent().first;
        if (extent.rank() != 1 || extent[0] != payload.size())
            throw std::runtime_error{ std::string{ __PRETTY_FUNCTION__ } + " - incompatible extent : gamma" };

        hdf5::Extent const start = { 0U, offsetof(Particle, gamma) / unit_size };
        hdf5::Extent const count = { payload.size(), sizeof(T) / unit_size };
        mspace.select(H5S_SELECT_SET, start, count);
        fspace.select_all();
        dset.read(fspace, payload.data(), type, mspace);
    }
    {
        auto const type   = hdf5::make_type<long>();
        using T           = std::decay_t<decltype(std::declval<Particle>().id)>;
        auto       dset   = root.dataset("id");
        auto       fspace = dset.space();
        auto const extent = fspace.simple_extent().first;
        if (extent.rank() != 1 || extent[0] != payload.size())
            throw std::runtime_error{ std::string{ __PRETTY_FUNCTION__ } + " - incompatible extent : id" };

        hdf5::Extent const start = { 0U, offsetof(Particle, id) / unit_size };
        hdf5::Extent const count = { payload.size(), sizeof(T) / unit_size };
        mspace.select(H5S_SELECT_SET, start, count);
        fspace.select_all();
        dset.read(fspace, payload.data(), type, mspace);
    }

    // distribute
    // TODO: The reverse operation can be removed.
    // This is to make the sequence the same as in the version of multi-thread particle loading.
    std::reverse(payload.begin(), payload.end());

    comm.bcast<long>(long(payload.size()), master).unpack([](auto) {});
    comm.bcast<Particle>(std::move(payload), master)
        .unpack(
            [](auto payload, PartSpecies &sp, bool append) {
                sp.load_ptls(std::move(payload), append);
            },
            sp, false);
}
long Snapshot::load_master(Domain &domain) const &
{
    // open hdf5 file and root group
    hdf5::Group root = hdf5::File(hdf5::File::rdonly_tag{}, filepath().c_str()).group("rel_pic_1d");

    // verify signature
    if (signature != root.attribute("signature").read<decltype(signature)>())
        throw std::runtime_error{ std::string{ __PRETTY_FUNCTION__ } + " - signature verification failed" };

    // B & E
    load_helper(root, domain.bfield, "bfield");
    load_helper(root, domain.efield, "efield");

    // particles
    for (unsigned i = 0; i < domain.part_species.size(); ++i) {
        PartSpecies      &sp    = domain.part_species.at(i);
        std::string const gname = std::string{ "part_species" } + '@' + std::to_string(i);
        auto const        group = root.group(gname.c_str());
        load_helper(group, sp);
    }

    // cold fluid
    for (unsigned i = 0; i < domain.cold_species.size(); ++i) {
        ColdSpecies      &sp    = domain.cold_species.at(i);
        std::string const gname = std::string{ "cold_species" } + '@' + std::to_string(i);
        auto const        group = root.group(gname.c_str());
        load_helper(group, sp.mom0_full, "mom0_full");
        load_helper(group, sp.mom1_full, "mom1_full");
    }

    // step count
    auto const step_count = root.attribute("step_count").read<long>();
    return comm.bcast<long>(step_count, master).unpack([step_count](auto) {
        // ignoring the arg is deliberate because bcasting to root itself does not work for some MPI impls
        return step_count;
    });
}
long Snapshot::load_worker(Domain &domain) const &
{
    // B & E
    comm.scatter<1>(nullptr, domain.bfield.begin(), domain.bfield.end(), master);
    comm.scatter<1>(nullptr, domain.efield.begin(), domain.efield.end(), master);

    // particles
    for (PartSpecies &sp : domain.part_species) {
        long const            Np = comm.bcast<long>({}, master);
        std::vector<Particle> payload(static_cast<unsigned long>(Np));
        comm.bcast<Particle>(payload, master)
            .unpack(
                [](auto payload, PartSpecies &sp, bool append) {
                    sp.load_ptls(std::move(payload), append);
                },
                sp, false);
    }

    // cold fluid
    for (ColdSpecies &sp : domain.cold_species) {
        comm.scatter<0>(nullptr, sp.mom0_full.begin(), sp.mom0_full.end(), master);
        comm.scatter<1>(nullptr, sp.mom1_full.begin(), sp.mom1_full.end(), master);
    }

    // step count
    return comm.bcast<long>({}, master);
}
PIC1D_END_NAMESPACE
