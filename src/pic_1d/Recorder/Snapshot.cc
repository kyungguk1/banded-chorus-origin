/*
 * Copyright (c) 2020-2021, Kyungguk Min
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * - Redistributions of source code must retain the above copyright notice, this
 *  list of conditions and the following disclaimer.
 *
 * - Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "Snapshot.h"

#include <algorithm>
#include <fstream>
#include <functional> // std::hash, std::mem_fn
#include <iterator>
#include <stdexcept>
#include <type_traits>

namespace {
template <class Tuple> struct Hash;
template <class T> Hash(T const &t) -> Hash<T>;
//
template <class... Ts> struct Hash<std::tuple<Ts...>> {
    std::tuple<Ts...> const t;

    [[nodiscard]] constexpr operator std::size_t() const noexcept { return operator()(); }
    [[nodiscard]] constexpr std::size_t operator()() const noexcept
    {
        return hash(std::index_sequence_for<Ts...>{});
    }

private:
    template <std::size_t... Is>
    [[nodiscard]] constexpr std::size_t hash(std::index_sequence<Is...>) const noexcept
    {
        std::size_t hash = 0;
        return (..., ((hash <<= 1) ^= this->hash(std::get<Is>(t)))), hash;
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
: comm{std::move(_comm), tag}, signature{Hash{serialize(params)}}, wd{params.working_directory}
{
    if (!comm->operator bool())
        throw std::invalid_argument{std::string{__PRETTY_FUNCTION__} + " - invalid mpi::Comm"};

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

template <class T, long N>
auto P1D::Snapshot::save_helper(hdf5::Group &root, GridQ<T, N> const &grid,
                                std::string const &basename) const -> hdf5::Dataset
{
    auto payload = *comm.gather<T>({grid.begin(), grid.end()}, master);

    // fixed dataset
    auto space = hdf5::Space::simple(payload.size());
    auto dset  = root.dataset(basename.c_str(), hdf5::make_type<T>(), space);

    // export
    space.select_all();
    dset.write(space, payload.data(), space);

    return dset;
}
auto P1D::Snapshot::save_helper(hdf5::Group &root, PartSpecies const &sp,
                                std::string const &basename) const -> hdf5::Dataset
{
    // chunked dataset
    auto dcpl = hdf5::PList::dcpl();
    dcpl.set_chunk(1024U);
    auto dset = root.dataset(basename.c_str(), hdf5::make_type<Particle>(),
                             hdf5::Space::simple(0U, nullptr), hdf5::PList::dapl(), dcpl);

    // export
    auto payload = sp.dump_ptls();
    auto tk      = comm.ibsend(std::move(payload), master);
    {
        unsigned long accum_count = 0;
        unsigned long start       = 0;
        for (int rank = 0, size = comm.size(); rank < size; ++rank, start = accum_count) {
            payload = comm.recv(std::move(payload), rank);
            accum_count += payload.size();

            auto mspace = hdf5::Space::simple(payload.size());
            mspace.select_all();

            dset.set_extent(accum_count);
            auto fspace = dset.space();
            fspace.select(H5S_SELECT_SET, start, payload.size());
            dset.write(fspace, payload.data(), mspace);
        }
    }
    std::move(tk).wait();

    return dset;
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
        PartSpecies const &sp     = domain.part_species.at(i);
        std::string const  prefix = std::string{"part_species_"} + std::to_string(i);
        save_helper(root, sp, prefix + "-particles") << sp;
    }

    // cold fluid
    for (unsigned i = 0; i < domain.cold_species.size(); ++i) {
        ColdSpecies const &sp     = domain.cold_species.at(i);
        std::string const  prefix = std::string{"cold_species_"} + std::to_string(i);
        save_helper(root, sp.mom0_full, prefix + "-mom0_full") << sp;
        save_helper(root, sp.mom1_full, prefix + "-mom1_full") << sp;
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

template <class T, long N>
auto P1D::Snapshot::load_helper(hdf5::Group const &root, GridQ<T, N> &grid,
                                std::string const &basename) const -> hdf5::Dataset
{
    std::vector<T> payload(static_cast<unsigned long>(grid.size() * comm.size()));

    // open dataset
    auto       dset   = root.dataset(basename.c_str());
    auto       space  = dset.space();
    auto const extent = space.simple_extent().first;
    if (extent.rank() != 1)
        throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                 + " - incompatible extent rank : " + basename};
    if (extent.front() != payload.size())
        throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                 + " - incompatible grid size : " + basename};

    // import
    space.select_all();
    dset.read(space, payload.data(), space);

    // distribute
    comm.scatter(payload.data(), grid.begin(), grid.end(), master);

    return dset;
}
auto P1D::Snapshot::load_helper(hdf5::Group const &root, PartSpecies &sp,
                                std::string const &basename) const -> hdf5::Dataset
{
    std::vector<Particle> payload;

    // open dataset
    auto       dset   = root.dataset(basename.c_str());
    auto       space  = dset.space();
    auto const extent = space.simple_extent().first;
    if (extent.rank() != 1)
        throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                 + " - incompatible extent rank : " + basename};
    if (extent.front() != static_cast<unsigned long>(sp->Nc * sp.params.Nx))
        throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                 + " - incompatible particle size : " + basename};
    payload.resize(extent.front());

    // import
    space.select_all();
    dset.read(space, payload.data(), space);

    // distribute
    std::reverse(payload.begin(), payload.end()); // This is to make the sequence the same as in
                                                  // the version of multi-thread particle loading.
    auto last = payload.crbegin();
    for (long i = 0; i < sp.params.Nx; ++i) { // assumption here is that the number of particles is
                                              // divisible to the number of grid points.
        std::advance(last, sp->Nc);
        comm.bcast<Particle>({payload.crbegin(), last}, master)
            .unpack(
                [](auto payload, PartSpecies &sp, bool append) {
                    sp.load_ptls(std::move(payload), append);
                },
                sp, i);
        payload.erase(last.base(), payload.end());
    }
    if (!payload.empty())
        throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                 + " - particles still remaining after distribution : " + basename};

    return dset;
}
long P1D::Snapshot::load_master(Domain &domain) const &
{
    // open hdf5 file and root group
    hdf5::Group root = hdf5::File(hdf5::File::rdonly_tag{}, filepath().c_str()).group("pic_1d");

    // verify signature
    if (signature != root.attribute("signature").read<decltype(signature)>())
        throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                 + " - signature verification failed"};

    // B & E
    load_helper(root, domain.bfield, "bfield");
    load_helper(root, domain.efield, "efield");

    // particles
    for (unsigned i = 0; i < domain.part_species.size(); ++i) {
        PartSpecies &     sp     = domain.part_species.at(i);
        std::string const prefix = std::string{"part_species_"} + std::to_string(i);
        load_helper(root, sp, prefix + "-particles");

        auto const count = *comm.all_reduce<long>(parallel::mpi::ReduceOp::plus(),
                                                  static_cast<long>(sp.bucket.size()));
        if (sp->Nc * sp.params.Nx != count)
            throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                     + " - particle count inconsistent for species "
                                     + std::to_string(i)};
    }

    // cold fluid
    for (unsigned i = 0; i < domain.cold_species.size(); ++i) {
        ColdSpecies &     sp     = domain.cold_species.at(i);
        std::string const prefix = std::string{"cold_species_"} + std::to_string(i);
        load_helper(root, sp.mom0_full, prefix + "-mom0_full");
        load_helper(root, sp.mom1_full, prefix + "-mom1_full");
    }

    // step count
    return comm.bcast<long>(root.attribute("step_count").read<long>(), master);
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
