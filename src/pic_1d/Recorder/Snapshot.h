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

#ifndef Snapshot_h
#define Snapshot_h

#include "../Core/Domain.h"
#include "../TypeMaps.h"

#include <ParallelKit/ParallelKit.h>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

PIC1D_BEGIN_NAMESPACE
namespace thread {
class Snapshot {
public:
    using message_dispatch_t
        = parallel::MessageDispatch<std::vector<Scalar>, std::vector<Vector>, std::vector<Tensor>,
                                    std::shared_ptr<std::vector<Particle> const>,
                                    std::vector<Particle>, long>;
    using interthread_comm_t = message_dispatch_t::Communicator;
    using ticket_t           = message_dispatch_t::Ticket;
    using rank_t             = message_dispatch_t::Rank;

    static message_dispatch_t dispatch;
    interthread_comm_t const  comm;
    unsigned const            size;
    static constexpr unsigned master = 0;
    [[nodiscard]] bool        is_master() const noexcept { return master == comm.rank; }

private:
    void (Snapshot::*save)(Domain const &domain) const &;
    long (Snapshot::*load)(Domain &domain) const &;
    long const          step_count;
    std::size_t const   signature;
    std::vector<rank_t> all_ranks;

    [[nodiscard]] static std::string filepath(std::string const &wd, std::string_view basename);

public:
    explicit Snapshot(unsigned rank, unsigned size, ParamSet const &params, long step_count);

private: // load/save
    void               save_master(Domain const &domain) const &;
    void               save_worker(Domain const &domain) const &;
    [[nodiscard]] long load_master(Domain &domain) const &;
    [[nodiscard]] long load_worker(Domain &domain) const &;

private: // load/save interface
    friend void operator<<(Snapshot &&snapshot, Domain const &domain)
    {
        return (snapshot.*snapshot.save)(domain);
    }
    [[nodiscard]] friend long operator>>(Snapshot &&snapshot, Domain &domain)
    {
        return (snapshot.*snapshot.load)(domain);
    }
};
} // namespace thread

namespace mpi {
class Snapshot {
public:
    using interprocess_comm_t = parallel::Communicator<Scalar, Vector, Tensor, Particle, long>;
    using rank_t              = parallel::mpi::Rank;

    static constexpr parallel::mpi::Tag tag{599};

private:
    interprocess_comm_t const comm;
    int                       size{};
    rank_t                    rank{-1};
    std::size_t const         signature;
    std::string const         wd; // working directory

    static constexpr rank_t   master{0};
    [[nodiscard]] bool        is_master() const noexcept { return master == this->rank; }
    [[nodiscard]] std::string filepath(std::string_view basename) const;

public:
    explicit Snapshot(parallel::mpi::Comm comm, ParamSet const &params);

private: // load/save
    void (Snapshot::*save)(Domain const &domain, long step_count) const &;
    long (Snapshot::*load)(Domain &domain) const &;

    template <class T, long N>
    void save_helper(GridQ<T, N> const &payload, long step_count, std::string_view basename) const;
    void save_helper(PartSpecies const &payload, long step_count, std::string_view basename) const;
    void save_master(Domain const &domain, long step_count) const &;
    void save_worker(Domain const &domain, long step_count) const &;

    template <class T, long N>
    long               load_helper(GridQ<T, N> &payload, std::string_view basename) const;
    long               load_helper(PartSpecies &payload, std::string_view basename) const;
    [[nodiscard]] long load_master(Domain &domain) const &;
    [[nodiscard]] long load_worker(Domain &domain) const &;

private: // load/save interface
    friend void save(Snapshot &&snapshot, Domain const &domain, long step_count)
    {
        return (snapshot.*snapshot.save)(domain, step_count);
    }
    [[nodiscard]] friend long load(Snapshot &&snapshot, Domain &domain)
    {
        return (snapshot.*snapshot.load)(domain);
    }
};
} // namespace mpi
PIC1D_END_NAMESPACE

#endif /* Snapshot_h */
