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

namespace {
template <class T, std::enable_if_t<std::is_trivially_copyable_v<T>, long> = 0L>
decltype(auto) write(std::ostream &os, T const &payload)
{
    return os.write(reinterpret_cast<char const *>(std::addressof(payload)), sizeof(T));
}
template <class T, std::enable_if_t<std::is_trivially_copyable_v<T>, long> = 0L>
decltype(auto) write(std::ostream &os, std::vector<T> const &payload)
{
    return os.write(reinterpret_cast<char const *>(payload.data()),
                    static_cast<long>(payload.size() * sizeof(T)));
}
//
template <class T, std::enable_if_t<std::is_trivially_copyable_v<T>, long> = 0L>
decltype(auto) read(std::istream &is, T &payload)
{
    return is.read(reinterpret_cast<char *>(std::addressof(payload)), sizeof(T));
}
template <class T, std::enable_if_t<std::is_trivially_copyable_v<T>, long> = 0L>
decltype(auto) read(std::istream &is, std::vector<T> &payload)
{
    return is.read(reinterpret_cast<char *>(payload.data()),
                   static_cast<long>(payload.size() * sizeof(T)));
}
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
std::string P1D::Snapshot::filepath(std::string_view const basename) const
{
    std::string const filename
        = std::string{"snapshot"} + "-" + std::string{basename} + ".snapshot";
    return wd + "/" + filename;
}

template <class T, long N>
void P1D::Snapshot::save_helper(GridQ<T, N> const &grid, long const step_count,
                                     std::string_view const basename) const
{
    std::string const path = filepath(basename);
    std::ofstream     os{path};
    if (!os)
        throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                 + " - file open failed : " + path};

    if (!write(os, signature))
        throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                 + " - writing signature failed : " + path};

    if (!write(os, step_count))
        throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                 + " - writing step count failed : " + path};

    if (!write(os, *comm.gather<T>({grid.begin(), grid.end()}, master)))
        throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                 + " - writing payload failed : " + path};
}
void P1D::Snapshot::save_helper(PartSpecies const &sp, long const step_count,
                                     std::string_view const basename) const
{
    std::string const path = filepath(basename);
    std::ofstream     os{path};
    if (!os)
        throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                 + " - file open failed : " + path};

    if (!write(os, signature))
        throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                 + " - writing signature failed : " + path};

    if (!write(os, step_count))
        throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                 + " - writing step count failed : " + path};

    auto payload = sp.dump_ptls();

    // particle count
    auto const count = *comm.all_reduce<long>(parallel::mpi::ReduceOp::plus(),
                                              static_cast<long>(payload.size()));
    if (!write(os, count))
        throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                 + " - writing particle count failed : " + path};

    // particle dump
    auto tk = comm.ibsend(std::move(payload), master);
    for (int rank = 0, size = comm.size(); rank < size; ++rank) {
        payload = comm.recv(std::move(payload), rank);
        if (!write(os, std::move(payload)))
            throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                     + " - writing particles failed : " + path + " with rank "
                                     + std::to_string(rank)};
    }
    std::move(tk).wait();
}
void P1D::Snapshot::save_master(Domain const &domain, long const step_count) const &
{
    // B & E
    save_helper(domain.bfield, step_count, "bfield");
    save_helper(domain.efield, step_count, "efield");

    // particles
    for (unsigned i = 0; i < domain.part_species.size(); ++i) {
        PartSpecies const &sp     = domain.part_species.at(i);
        std::string const  prefix = std::string{"part_species_"} + std::to_string(i);
        save_helper(sp, step_count, prefix + "-particles");
    }

    // cold fluid
    for (unsigned i = 0; i < domain.cold_species.size(); ++i) {
        ColdSpecies const &sp     = domain.cold_species.at(i);
        std::string const  prefix = std::string{"cold_species_"} + std::to_string(i);
        save_helper(sp.mom0_full, step_count, prefix + "-mom0_full");
        save_helper(sp.mom1_full, step_count, prefix + "-mom1_full");
    }
}
void P1D::Snapshot::save_worker(Domain const &domain, long) const &
{
    // B & E
    comm.gather<1>(domain.bfield.begin(), domain.bfield.end(), nullptr, master);
    comm.gather<1>(domain.efield.begin(), domain.efield.end(), nullptr, master);

    // particles
    for (PartSpecies const &sp : domain.part_species) {
        auto payload = sp.dump_ptls();
        comm.all_reduce<long>(parallel::mpi::ReduceOp::plus(), static_cast<long>(payload.size()))
            .unpack([](auto) {});
        comm.ibsend(std::move(payload), master).wait();
    }

    // cold fluid
    for (ColdSpecies const &sp : domain.cold_species) {
        comm.gather<0>(sp.mom0_full.begin(), sp.mom0_full.end(), nullptr, master);
        comm.gather<1>(sp.mom1_full.begin(), sp.mom1_full.end(), nullptr, master);
    }
}

template <class T, long N>
long P1D::Snapshot::load_helper(GridQ<T, N> &grid, std::string_view const basename) const
{
    std::string const path = filepath(basename);
    std::ifstream     is{path};
    if (!is)
        throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                 + " - file open failed : " + path};

    std::size_t signature;
    if (!read(is, signature))
        throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                 + " - reading signature failed : " + path};
    if (this->signature != signature)
        throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                 + " - incompatible signature : " + path};

    long step_count;
    if (!read(is, step_count))
        throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                 + " - reading step count failed : " + path};

    std::vector<T> payload(static_cast<unsigned long>(grid.size() * comm.size()));
    if (!read(is, payload))
        throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                 + " - reading payload failed : " + path};
    if (char dummy; !read(is, dummy).eof())
        throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                 + " - payload not fully read : " + path};

    // distribute payload
    comm.scatter(payload.data(), grid.begin(), grid.end(), master);

    return step_count;
}
long P1D::Snapshot::load_helper(PartSpecies &sp, std::string_view const basename) const
{
    std::string const path = filepath(basename);
    std::ifstream     is{path};
    if (!is)
        throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                 + " - file open failed : " + path};

    std::size_t signature;
    if (!read(is, signature))
        throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                 + " - reading signature failed : " + path};
    if (this->signature != signature)
        throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                 + " - incompatible signature : " + path};

    long step_count;
    if (!read(is, step_count))
        throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                 + " - reading step count failed : " + path};

    // particle count
    long ptl_count;
    if (!read(is, ptl_count))
        throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                 + " - reading particle count failed : " + path};
    if (sp->Nc * sp.params.Nx != ptl_count)
        throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                 + " - particle count read inconsistent : " + path};

    // particle load
    std::vector<Particle> payload(static_cast<unsigned long>(ptl_count));
    if (!read(is, payload))
        throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                 + " - reading particles failed : " + path};
    if (char dummy; !read(is, dummy).eof())
        throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                 + " - particles not fully read : " + path};

    // This is to make the sequence the same as in the version of multi-thread particle loading.
    std::reverse(payload.begin(), payload.end());

    // particle distribute
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
                                 + " - particles still remaining after distribution : " + path};

    return step_count;
}
long P1D::Snapshot::load_master(Domain &domain) const &
{
    // B & E
    long const step_count = load_helper(domain.bfield, "bfield");
    if (step_count != load_helper(domain.efield, "efield"))
        throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                 + " - incompatible step_count : efield"};

    // particles
    for (unsigned i = 0; i < domain.part_species.size(); ++i) {
        PartSpecies &     sp     = domain.part_species.at(i);
        std::string const prefix = std::string{"part_species_"} + std::to_string(i);
        if (step_count != load_helper(sp, prefix + "-particles"))
            throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                     + " - incompatible step_count : " + std::to_string(i)
                                     + "th part_species"};

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
        if (step_count != load_helper(sp.mom0_full, prefix + "-mom0_full"))
            throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                     + " - incompatible step_count : " + std::to_string(i)
                                     + "th cold_mom<0>"};
        if (step_count != load_helper(sp.mom1_full, prefix + "-mom1_full"))
            throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                     + " - incompatible step_count : " + std::to_string(i)
                                     + "th cold_mom<1>"};
    }

    // step count
    return comm.bcast(step_count, master);
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
