//
// Copyright (c) 2019, Kyungguk Min
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// - Redistributions of source code must retain the above copyright notice, this
//  list of conditions and the following disclaimer.
//
// - Redistributions in binary form must reproduce the above copyright notice,
//  this list of conditions and the following disclaimer in the documentation
//  and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#ifndef MasterDelegate_h
#define MasterDelegate_h

#include "WorkerDelegate.h"

#include <array>
#include <type_traits>

HYBRID1D_BEGIN_NAMESPACE
class MasterDelegate final : public Delegate {
    using ticket_t = WorkerDelegate::message_dispatch_t::Ticket;

public:
    std::array<WorkerDelegate, ParamSet::number_of_particle_parallelism - 1> workers{};
    mutable // access of methods in message dispatcher is thread-safe
        WorkerDelegate::message_dispatch_t dispatch{
            ParamSet::number_of_particle_parallelism}; // each master thread in domain decomposition
                                                       // must have its own message dispatcher
    WorkerDelegate::interthread_comm_t comm{};
    Delegate *const                    delegate; // serial version
    std::vector<unsigned>              all_but_master;

    ~MasterDelegate() override;
    explicit MasterDelegate(Delegate *delegate);

private:
    void once(Domain &) const override;
    void prologue(Domain const &, long) const override;
    void epilogue(Domain const &, long) const override;
    void pass(Domain const &, PartSpecies &) const override;
    void pass(Domain const &, ColdSpecies &) const override;
    void pass(Domain const &, BField &) const override;
    void pass(Domain const &, EField &) const override;
    void pass(Domain const &, Charge &) const override;
    void pass(Domain const &, Current &) const override;
    void gather(Domain const &, Charge &) const override;
    void gather(Domain const &, Current &) const override;
    void gather(Domain const &, PartSpecies &) const override;

private: // helpers
    template <class T, long N> void broadcast_to_workers(GridQ<T, N> const &payload) const;
    template <class T, long N> void collect_from_workers(GridQ<T, N> &buffer) const;

public: // wrap the loop with setup/teardown logic included
    template <class F, class... Args> [[nodiscard]] auto wrap_loop(F &&f, Args &&...args)
    {
        return [this, f, args...](Domain *domain) mutable { // intentional capture by copy
            setup(*domain);
            std::invoke(std::forward<F>(f), std::move(args)...); // hence move is used
            teardown(*domain);
        };
    }

private:
    void collect(Domain const &, PartSpecies &) const;
    void distribute(Domain const &, PartSpecies &) const;

public:
    void setup(Domain &) const;
    void teardown(Domain &) const;
};
HYBRID1D_END_NAMESPACE

#endif /* MasterDelegate_h */
