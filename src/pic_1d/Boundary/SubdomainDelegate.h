/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include "Delegate.h"
#include <PIC/TypeMaps.h>

#include <ParallelKit/ParallelKit.h>

PIC1D_BEGIN_NAMESPACE
class SubdomainDelegate : public Delegate {
    using interprocess_comm_t = parallel::Communicator<Scalar, Vector, Tensor, Particle>;
    using rank_t              = parallel::mpi::Rank;

    interprocess_comm_t comm;
    rank_t              left_{ -1 };
    rank_t              right{ -1 };

    static constexpr rank_t master{ 0 };
    [[nodiscard]] bool      is_master() const { return master == comm->rank(); }

    // these must be consistent with the definition in ParamSet
    [[nodiscard]] bool is_leftmost_subdomain() const { return comm->rank() == 0; }
    [[nodiscard]] bool is_rightmost_subdomain() const { return comm->rank() == comm.size() - 1; }

public:
    explicit SubdomainDelegate(parallel::mpi::Comm comm);

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
    template <class T, long N>
    void pass(Grid<T, N, Pad> &) const;
    template <class T, long N>
    void gather(Grid<T, N, Pad> &) const;
};
PIC1D_END_NAMESPACE
