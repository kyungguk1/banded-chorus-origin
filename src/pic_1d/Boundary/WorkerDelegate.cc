/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "WorkerDelegate.h"
#include "MasterDelegate.h"

PIC1D_BEGIN_NAMESPACE
void WorkerDelegate::setup(Domain &domain) const
{
    // distribute particles to workers
    //
    for (PartSpecies &sp : domain.part_species) {
        sp.Nc /= ParamSet::number_of_particle_parallelism;
        distribute(domain, sp);
    }
}
void WorkerDelegate::distribute(Domain const &, PartSpecies &sp) const
{
    // distribute particles to workers
    //
    sp.bucket = comm.recv<decltype(sp.bucket)>(master->comm.rank);
}

void WorkerDelegate::teardown(Domain &domain) const
{
    // collect particles to master
    //
    for (PartSpecies &sp : domain.part_species) {
        collect(domain, sp);
        sp.Nc *= ParamSet::number_of_particle_parallelism;
    }
}
void WorkerDelegate::collect(Domain const &, PartSpecies &sp) const
{
    // collect particles to master
    //
    comm.send(std::move(sp.bucket), master->comm.rank).wait();
}

void WorkerDelegate::prologue(Domain const &domain, long const i) const
{
    master->delegate->prologue(domain, i);
}
void WorkerDelegate::epilogue(Domain const &domain, long const i) const
{
    master->delegate->epilogue(domain, i);
}
void WorkerDelegate::once(Domain &domain) const
{
    master->delegate->once(domain);
}
void WorkerDelegate::pass(Domain const &, PartSpecies &sp)
{
    auto &[L, R] = buckets.cleared(); // be careful not to access it from multiple threads
                                      // be sure to clear the contents before use
    master->delegate->partition(sp, L, R);
    //
    comm.recv<0>(master->comm.rank).unpack([&L = L, &R = R](auto payload) {
        payload.first->swap(L);
        payload.second->swap(R);
    });
    //
    sp.bucket.insert(sp.bucket.cend(), L.cbegin(), L.cend());
    sp.bucket.insert(sp.bucket.cend(), R.cbegin(), R.cend());
}
void WorkerDelegate::pass(Domain const &, ColdSpecies &sp) const
{
    recv_from_master(sp.mom0_full);
    recv_from_master(sp.mom1_full);
}
void WorkerDelegate::pass(Domain const &, BField &bfield) const
{
    recv_from_master(bfield);
}
void WorkerDelegate::pass(Domain const &, EField &efield) const
{
    recv_from_master(efield);
}
void WorkerDelegate::pass(Domain const &, Current &current) const
{
    recv_from_master(current);
}
void WorkerDelegate::gather(Domain const &, Current &current) const
{
    reduce_to_master(current);
    recv_from_master(current);
}
void WorkerDelegate::gather(Domain const &, PartSpecies &sp) const
{
    {
        reduce_to_master(sp.moment<0>());
        reduce_to_master(sp.moment<1>());
        reduce_to_master(sp.moment<2>());
    }
    {
        recv_from_master(sp.moment<0>());
        recv_from_master(sp.moment<1>());
        recv_from_master(sp.moment<2>());
    }
}

template <class T, long N>
void WorkerDelegate::recv_from_master(Grid<T, N, Pad> &buffer) const
{
    comm.recv<Grid<T, N, Pad> const *>(master->comm.rank).unpack([&buffer](auto payload) {
        buffer = *payload;
    });
}
template <class T, long N>
void WorkerDelegate::reduce_to_master(Grid<T, N, Pad> const &payload) const
{
    comm.send(&payload, master->comm.rank).wait(); // must wait for delivery receipt
}
PIC1D_END_NAMESPACE
