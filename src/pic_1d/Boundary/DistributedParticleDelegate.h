/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include "Delegate.h"
#include <PIC/TypeMaps.h>

#include <ParallelKit/ParallelKit.h>
#include <array>

PIC1D_BEGIN_NAMESPACE
class DistributedParticleDelegate : public Delegate {
    using interprocess_comm_t = parallel::Communicator<Scalar, Vector, Tensor>;
    using rank_t              = parallel::mpi::Rank;

    interprocess_comm_t                    comm;
    Delegate const                        *subdomain_delegate;
    std::array<parallel::mpi::ReduceOp, 3> reduce_plus;

public:
    DistributedParticleDelegate(parallel::mpi::Comm comm, Delegate const *subdomain_delegate);

private:
    void once(Domain &) const override;
    void prologue(Domain const &, long) const override;
    void epilogue(Domain const &, long) const override;
    void pass(Domain const &, PartSpecies &) const override;
    void pass(Domain const &, ColdSpecies &) const override;
    void pass(Domain const &, BField &) const override;
    void pass(Domain const &, EField &) const override;
    void pass(Domain const &, Current &) const override;
    void gather(Domain const &, Current &) const override;
    void gather(Domain const &, Species &) const override;

    // helpers
    template <unsigned I, class T, long S>
    void accumulate_distribute(Grid<T, S, Pad> &grid) const;
};
PIC1D_END_NAMESPACE
