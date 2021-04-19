/*
 * Copyright (c) 2019-2021, Kyungguk Min
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

#ifndef SubdomainDelegate_h
#define SubdomainDelegate_h

#include "./Delegate.h"
#include "./TypeMaps.h"

#include <ParallelKit/ParallelKit.h>

PIC1D_BEGIN_NAMESPACE
namespace thread {
class SubdomainDelegate : public Delegate {
public:
    using message_dispatch_t
        = parallel::MessageDispatch<Scalar const *, Vector const *, Tensor const *, PartBucket>;
    using interthread_comm_t = message_dispatch_t::Communicator;

    static message_dispatch_t dispatch;
    interthread_comm_t const  comm;
    unsigned const            size;
    unsigned const            left_;
    unsigned const            right;
    static constexpr unsigned master = 0;
    [[nodiscard]] bool        is_master() const noexcept { return master == comm.rank; }

public:
    SubdomainDelegate(unsigned rank, unsigned size);

private:
    void once(Domain &) const override;
    void prologue(Domain const &, long) const override {}
    void epilogue(Domain const &, long) const override {}

    // default implementation is periodic boundary condition
    //
    void pass(Domain const &, PartBucket &L_bucket, PartBucket &R_bucket) const override;
    void pass(Domain const &, ColdSpecies &) const override;
    void pass(Domain const &, BField &) const override;
    void pass(Domain const &, EField &) const override;
    void pass(Domain const &, Current &) const override;
    void gather(Domain const &, Current &) const override;
    void gather(Domain const &, PartSpecies &) const override;

private: // helpers
    template <class T, long N> void pass(GridQ<T, N> &) const;
    template <class T, long N> void gather(GridQ<T, N> &) const;
};
} // namespace thread

namespace mpi {
class SubdomainDelegate : public Delegate {
public:
    using interprocess_comm_t = parallel::Communicator<Scalar, Vector, Tensor, Particle>;
    using rank_t              = parallel::mpi::Rank;

    static constexpr parallel::mpi::Tag tag{598};

private:
    interprocess_comm_t comm;
    int                 size{};
    rank_t              rank{-1};
    rank_t              left_{-1};
    rank_t              right{-1};

    static constexpr rank_t master{0};
    [[nodiscard]] bool      is_master() const noexcept { return master == this->rank; }

public:
    SubdomainDelegate(parallel::mpi::Comm comm);

private:
    void once(Domain &) const override;
    void prologue(Domain const &, long) const override {}
    void epilogue(Domain const &, long) const override {}

    // default implementation is periodic boundary condition
    //
    void pass(Domain const &, PartBucket &L_bucket, PartBucket &R_bucket) const override;
    void pass(Domain const &, ColdSpecies &) const override;
    void pass(Domain const &, BField &) const override;
    void pass(Domain const &, EField &) const override;
    void pass(Domain const &, Current &) const override;
    void gather(Domain const &, Current &) const override;
    void gather(Domain const &, PartSpecies &) const override;

private: // helpers
    template <class T, long N> void pass(GridQ<T, N> &) const;
    template <class T, long N> void gather(GridQ<T, N> &) const;
};
} // namespace mpi
PIC1D_END_NAMESPACE

#endif /* SubdomainDelegate_h */
