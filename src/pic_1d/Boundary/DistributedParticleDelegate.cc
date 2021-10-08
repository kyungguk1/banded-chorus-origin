/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "DistributedParticleDelegate.h"

#include <stdexcept>
#include <utility>

PIC1D_BEGIN_NAMESPACE
DistributedParticleDelegate::DistributedParticleDelegate(parallel::mpi::Comm _comm, Delegate const *subdomain_delegate)
: comm{ std::move(_comm) }, subdomain_delegate{ subdomain_delegate }
{
    if (!comm->operator bool())
        throw std::invalid_argument{ __PRETTY_FUNCTION__ };

    using parallel::mpi::ReduceOp;
    reduce_plus = {
        ReduceOp::plus<Scalar>(true),
        ReduceOp::plus<Vector>(true),
        ReduceOp::plus<Tensor>(true),
    };
}

void DistributedParticleDelegate::prologue(Domain const &domain, long const i) const
{
    subdomain_delegate->prologue(domain, i);
}
void DistributedParticleDelegate::epilogue(Domain const &domain, long const i) const
{
    subdomain_delegate->epilogue(domain, i);
}
void DistributedParticleDelegate::once(Domain &domain) const
{
    subdomain_delegate->once(domain);
}
void DistributedParticleDelegate::pass(Domain const &domain, PartSpecies &sp) const
{
    subdomain_delegate->pass(domain, sp);
}
void DistributedParticleDelegate::pass(Domain const &domain, ColdSpecies &sp) const
{
    subdomain_delegate->pass(domain, sp);
}
void DistributedParticleDelegate::pass(Domain const &domain, BField &bfield) const
{
    subdomain_delegate->pass(domain, bfield);
}
void DistributedParticleDelegate::pass(Domain const &domain, EField &efield) const
{
    subdomain_delegate->pass(domain, efield);
}
void DistributedParticleDelegate::pass(Domain const &domain, Current &current) const
{
    subdomain_delegate->pass(domain, current);
}
void DistributedParticleDelegate::gather(Domain const &domain, Current &current) const
{
    subdomain_delegate->gather(domain, current);
    accumulate_distribute<1>(current);
}
void DistributedParticleDelegate::gather(Domain const &domain, Species &sp) const
{
    subdomain_delegate->gather(domain, sp);
    {
        accumulate_distribute<0>(sp.moment<0>());
        accumulate_distribute<1>(sp.moment<1>());
        accumulate_distribute<2>(sp.moment<2>());
    }
}

template <unsigned I, class T, long S>
void DistributedParticleDelegate::accumulate_distribute(Grid<T, S, Pad> &grid) const
{
    comm.all_reduce<I>(std::get<I>(reduce_plus), grid.begin(), grid.end());
}
PIC1D_END_NAMESPACE
