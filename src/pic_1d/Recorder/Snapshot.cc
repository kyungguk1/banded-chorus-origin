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

// MARK:- thread::Snapshot
//
P1D::thread::Snapshot::message_dispatch_t P1D::thread::Snapshot::dispatch{
    P1D::ParamSet::number_of_subdomains};
P1D::thread::Snapshot::Snapshot(unsigned const rank, unsigned const size, ParamSet const &params,
                                long const step_count)
: comm{dispatch.comm(rank)}
, size{size}
, step_count{step_count}
, signature{Hash{serialize(params)}}
, all_ranks{}
{
    if (size > ParamSet::number_of_subdomains)
        throw std::invalid_argument{__PRETTY_FUNCTION__};

    // method dispatch
    //
    if (is_master()) {
        save = &Snapshot::save_master;
        load = &Snapshot::load_master;
    } else {
        save = &Snapshot::save_worker;
        load = &Snapshot::load_worker;
    }

    // participants
    //
    for (unsigned rank = 0; is_master() && rank < size; ++rank) {
        all_ranks.emplace_back(rank);
    }
}

std::string P1D::thread::Snapshot::filepath(std::string const &wd, std::string_view const basename)
{
    std::string const filename
        = std::string{"snapshot"} + "-" + std::string{basename} + ".snapshot";
    return wd + "/" + filename;
}

// MARK:- Save
//
namespace {
template <class T, long N> [[nodiscard]] std::vector<T> pack(P1D::GridQ<T, N> const &payload)
{
    return {payload.begin(), payload.end()};
}
[[nodiscard]] auto pack(P1D::PartSpecies const &sp)
{
    return sp.dump_ptls();
}
//
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
} // namespace
void P1D::thread::Snapshot::save_master(Domain const &domain) const &
{
    auto save_grid = [this, wd = domain.params.working_directory](auto const &           payload,
                                                                  std::string_view const basename) {
        std::string const path = filepath(wd, basename);
        if (std::ofstream os{path}; os) {
            if (!write(os, signature)) {
                throw std::runtime_error{path + " - writing signature failed"};
            }
            if (!write(os, step_count)) {
                throw std::runtime_error{path + " - writing step count failed"};
            }
            //
            auto tk = comm.send(pack(payload), master);
            comm.for_each<decltype(pack(payload))>(all_ranks, [&os, path, basename](auto payload) {
                if (!write(os, std::move(payload))) {
                    throw std::runtime_error{
                        path + " - writing payload failed : " + std::string{basename}};
                }
            });
            std::move(tk).wait();
        } else {
            throw std::runtime_error{path + " - file open failed"};
        }
    };
    auto save_ptls = [this, wd = domain.params.working_directory](PartSpecies const &    sp,
                                                                  std::string_view const basename) {
        std::string const path = filepath(wd, basename);
        if (std::ofstream os{path}; os) {
            if (!write(os, signature)) {
                throw std::runtime_error{path + " - writing signature failed"};
            }
            if (!write(os, step_count)) {
                throw std::runtime_error{path + " - writing step count failed"};
            }
            auto payload = pack(sp);
            // particle count
            auto tk1 = comm.send(static_cast<long>(payload.size()), master);
            if (!write(os, comm.reduce<long>(all_ranks, long{}, std::plus{}))) {
                throw std::runtime_error{
                    path + " - writing particle count failed : " + std::string{basename}};
            }
            std::move(tk1).wait();
            // particle dump
            auto tk2 = comm.send(std::move(payload), master);
            comm.for_each<decltype(payload)>(all_ranks, [&os, path, basename](auto payload) {
                if (!write(os, std::move(payload))) {
                    throw std::runtime_error{
                        path + " - writing particles failed : " + std::string{basename}};
                }
            });
            std::move(tk2).wait();
        } else {
            throw std::runtime_error{path + " - file open failed"};
        }
    };

    // B & E
    save_grid(domain.bfield, "bfield");
    save_grid(domain.efield, "efield");

    // particles
    for (unsigned i = 0; i < domain.part_species.size(); ++i) {
        PartSpecies const &sp     = domain.part_species.at(i);
        std::string const  prefix = std::string{"part_species_"} + std::to_string(i);
        save_ptls(sp, prefix + "-particles");
    }

    // cold fluid
    for (unsigned i = 0; i < domain.cold_species.size(); ++i) {
        ColdSpecies const &sp     = domain.cold_species.at(i);
        std::string const  prefix = std::string{"cold_species_"} + std::to_string(i);
        save_grid(sp.mom0_full, prefix + "-mom0_full");
        save_grid(sp.mom1_full, prefix + "-mom1_full");
    }
}
void P1D::thread::Snapshot::save_worker(
    Domain const &domain) const & // just wait because not a performace critical section
{
    // B & E
    comm.send(pack(domain.bfield), master).wait();
    comm.send(pack(domain.efield), master).wait();

    // particles
    for (PartSpecies const &sp : domain.part_species) {
        auto payload = pack(sp);
        comm.send(static_cast<long>(payload.size()), master).wait();
        comm.send(std::move(payload), master).wait(); // potential memory exhaustion if not wait
    }

    // cold fluid
    for (ColdSpecies const &sp : domain.cold_species) {
        comm.send(pack(sp.mom0_full), master).wait();
        comm.send(pack(sp.mom1_full), master).wait();
    }
}

// MARK: Load
//
namespace {
template <class T, long N>
void unpack_grid(std::vector<T>    payload,
                 P1D::GridQ<T, N> &to) noexcept(std::is_nothrow_move_assignable_v<T>)
{
    std::move(begin(payload), end(payload), to.begin());
}
template <class T> void unpack_ptls(std::shared_ptr<T const> payload, P1D::PartSpecies &sp)
{
    sp.load_ptls(*payload);
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
long P1D::thread::Snapshot::load_master(Domain &domain) const &
{
    long step_count{};
    auto load_grid = [this, wd = domain.params.working_directory,
                      &step_count](auto &to, std::string_view const basename) {
        std::string const path = filepath(wd, basename);
        if (std::ifstream is{path}; is) {
            if (std::size_t signature; !read(is, signature)) {
                throw std::runtime_error{path + " - reading signature failed"};
            } else if (this->signature != signature) {
                throw std::runtime_error{path + " - incompatible signature"};
            }
            if (!read(is, step_count)) {
                throw std::runtime_error{path + " - reading step count failed"};
            }
            //
            std::vector<decltype(pack(to))> payloads;
            payloads.reserve(all_ranks.size());
            for ([[maybe_unused]] rank_t const &rank : all_ranks) {
                if (!read(is, payloads.emplace_back(to.size()))) {
                    throw std::runtime_error{
                        path + " - reading payload failed : " + std::string{basename}};
                }
            }
            if (char dummy; !read(is, dummy).eof()) {
                throw std::runtime_error{path + " - payload not fully read"};
            }
            //
            auto tks = comm.scatter(std::move(payloads), all_ranks);
            unpack_grid(*comm.recv<decltype(pack(to))>(master), to);
            std::for_each(std::make_move_iterator(begin(tks)), std::make_move_iterator(end(tks)),
                          std::mem_fn(&ticket_t::wait));
        } else {
            throw std::runtime_error{path + " - file open failed"};
        }
    };
    auto load_ptls = [this, wd = domain.params.working_directory,
                      &step_count](PartSpecies &sp, std::string_view const basename) {
        std::string const path = filepath(wd, basename);
        if (std::ifstream is{path}; is) {
            if (std::size_t signature; !read(is, signature)) {
                throw std::runtime_error{path + " - reading signature failed"};
            } else if (this->signature != signature) {
                throw std::runtime_error{path + " - incompatible signature"};
            }
            if (!read(is, step_count)) {
                throw std::runtime_error{path + " - reading step count failed"};
            }
            std::shared_ptr<decltype(pack(sp))> payload;
            // particle count
            if (long size{}; !read(is, size)) {
                throw std::runtime_error{path + " - reading particle count failed"};
            } else if (sp->Nc * sp.params.Nx != size) {
                throw std::runtime_error{path + " - particle count read inconsistent"};
            } else {
                payload = std::make_shared<decltype(pack(sp))>(size);
            }
            // particle load
            if (!payload) {
                throw std::runtime_error{path + " - particle bucket not initialized"};
            } else if (!read(is, *payload)) {
                throw std::runtime_error{path + " - reading particles failed"};
            } else if (char dummy; !read(is, dummy).eof()) {
                throw std::runtime_error{path + " - particles not fully read"};
            }
            // sent payload must be alive until all workers got their particles loaded
            auto tks = comm.bcast<3>(payload, all_ranks);
            unpack_ptls(*comm.recv<3>(master), sp);
            std::for_each(std::make_move_iterator(begin(tks)), std::make_move_iterator(end(tks)),
                          std::mem_fn(&ticket_t::wait));
        } else {
            throw std::runtime_error{path + " - file open failed"};
        }
    };

    // B & E
    load_grid(domain.bfield, "bfield");
    load_grid(domain.efield, "efield");

    // particles
    for (unsigned i = 0; i < domain.part_species.size(); ++i) {
        PartSpecies &     sp     = domain.part_species.at(i);
        std::string const prefix = std::string{"part_species_"} + std::to_string(i);
        load_ptls(sp, prefix + "-particles");
        //
        auto tk = comm.send(static_cast<long>(sp.bucket.size()), master);
        if (sp->Nc * sp.params.Nx != comm.reduce<long>(all_ranks, long{}, std::plus{})) {
            throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                     + " - particle count inconsistent for species "
                                     + std::to_string(i)};
        }
        std::move(tk).wait();
    }

    // cold fluid
    for (unsigned i = 0; i < domain.cold_species.size(); ++i) {
        ColdSpecies &     sp     = domain.cold_species.at(i);
        std::string const prefix = std::string{"cold_species_"} + std::to_string(i);
        load_grid(sp.mom0_full, prefix + "-mom0_full");
        load_grid(sp.mom1_full, prefix + "-mom1_full");
    }

    // step count
    auto tks = comm.bcast(step_count, all_ranks);
    return comm.recv<long>(master);
    // std::for_each(std::make_move_iterator(begin(tks)), std::make_move_iterator(end(tks)),
    // std::mem_fn(&ticket_t::wait));
}
long P1D::thread::Snapshot::load_worker(Domain &domain) const &
{
    // B & E
    unpack_grid(*comm.recv<1>(master), domain.bfield);
    unpack_grid(*comm.recv<1>(master), domain.efield);

    // particles
    for (PartSpecies &sp : domain.part_species) {
        // received payload must be alive until all workers got their particles loaded
        unpack_ptls(*comm.recv<3>(master), sp);
        //
        comm.send(static_cast<long>(sp.bucket.size()), master).wait();
    }

    // cold fluid
    for (ColdSpecies &sp : domain.cold_species) {
        unpack_grid(*comm.recv<0>(master), sp.mom0_full);
        unpack_grid(*comm.recv<1>(master), sp.mom1_full);
    }

    // step count
    return comm.recv<long>(master);
}

// MARK:- mpi::Snapshot
//
P1D::mpi::Snapshot::Snapshot(parallel::mpi::Comm _comm, ParamSet const &params)
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
std::string P1D::mpi::Snapshot::filepath(std::string_view const basename) const
{
    std::string const filename
        = std::string{"snapshot"} + "-" + std::string{basename} + ".snapshot";
    return wd + "/" + filename;
}

template <class T, long N>
void P1D::mpi::Snapshot::save_helper(GridQ<T, N> const &grid, long const step_count,
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
void P1D::mpi::Snapshot::save_helper(PartSpecies const &sp, long const step_count,
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
void P1D::mpi::Snapshot::save_master(Domain const &domain, long const step_count) const &
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
void P1D::mpi::Snapshot::save_worker(Domain const &domain, long) const &
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
long P1D::mpi::Snapshot::load_helper(GridQ<T, N> &grid, std::string_view const basename) const
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
long P1D::mpi::Snapshot::load_helper(PartSpecies &sp, std::string_view const basename) const
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
long P1D::mpi::Snapshot::load_master(Domain &domain) const &
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
long P1D::mpi::Snapshot::load_worker(Domain &domain) const &
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
