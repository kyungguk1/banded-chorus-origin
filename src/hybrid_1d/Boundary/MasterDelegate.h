/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include "WorkerDelegate.h"

#include <array>
#include <type_traits>

HYBRID1D_BEGIN_NAMESPACE
class MasterDelegate final : public Delegate {
    using ticket_t = WorkerDelegate::message_dispatch_t::Ticket;
    using rank_t   = WorkerDelegate::message_dispatch_t::Rank;

public:
    std::array<WorkerDelegate, ParamSet::number_of_particle_parallelism - 1> workers{};
    mutable // access of methods in message dispatcher is thread-safe
        WorkerDelegate::message_dispatch_t dispatch{
            ParamSet::number_of_particle_parallelism
        }; // each master thread in domain decomposition must have its own message dispatcher
    WorkerDelegate::interthread_comm_t comm{};
    Delegate *const                    delegate; // serial version
    std::vector<rank_t>                all_but_master;

    ~MasterDelegate() override;
    explicit MasterDelegate(Delegate *delegate);

private:
    void once(Domain &) const override;
    void prologue(Domain const &, long) const override;
    void epilogue(Domain const &, long) const override;
    void pass(Domain const &, PartSpecies &) override;
    void pass(Domain const &, ColdSpecies &) const override;
    void pass(Domain const &, BField &) const override;
    void pass(Domain const &, EField &) const override;
    void pass(Domain const &, Charge &) const override;
    void pass(Domain const &, Current &) const override;
    void gather(Domain const &, Charge &) const override;
    void gather(Domain const &, Current &) const override;
    void gather(Domain const &, Species &) const override;

    // helpers
    template <class T, long N>
    void broadcast_to_workers(Grid<T, N, Pad> const &payload) const;
    template <class T, long N>
    void collect_from_workers(Grid<T, N, Pad> &buffer) const;

public: // wrap the loop with setup/teardown logic included
    template <class F, class... Args>
    [[nodiscard]] auto wrap_loop(F &&f, Args &&...args)
    {
        return [this, f, args...](Domain *domain) mutable { // intentional capture by copy
            setup(*domain);
            std::invoke(std::forward<F>(f), std::move(args)...); // hence, move is used
            teardown(*domain);
        };
    }

    void setup(Domain &) const;
    void teardown(Domain &) const;

private:
    void collect(Domain const &, PartSpecies &) const;
    void collect(Domain const &, ColdSpecies &) const;
    void distribute(Domain const &, PartSpecies &) const;
    void distribute(Domain const &, ColdSpecies &) const;
};
HYBRID1D_END_NAMESPACE
