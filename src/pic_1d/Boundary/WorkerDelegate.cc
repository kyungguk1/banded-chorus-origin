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

#include "WorkerDelegate.h"

#include "./MasterDelegate.h"

void P1D::WorkerDelegate::setup(Domain &domain) const
{
    // distribute particles to workers
    //
    for (PartSpecies &sp : domain.part_species) {
        sp.Nc /= ParamSet::number_of_particle_parallelism;
        distribute(domain, sp);
    }
}
void P1D::WorkerDelegate::distribute(Domain const &, PartSpecies &sp) const
{
    // distribute particles to workers
    //
    sp.bucket = comm.recv<PartBucket>(master->comm.rank);
}

void P1D::WorkerDelegate::teardown(Domain &domain) const
{
    // collect particles to master
    //
    for (PartSpecies &sp : domain.part_species) {
        collect(domain, sp);
        sp.Nc *= ParamSet::number_of_particle_parallelism;
    }
}
void P1D::WorkerDelegate::collect(Domain const &, PartSpecies &sp) const
{
    // collect particles to master
    //
    comm.send(std::move(sp.bucket), master->comm.rank).wait();
}

void P1D::WorkerDelegate::prologue(Domain const &domain, long const i) const
{
    master->delegate->prologue(domain, i);
}
void P1D::WorkerDelegate::epilogue(Domain const &domain, long const i) const
{
    master->delegate->epilogue(domain, i);
}
void P1D::WorkerDelegate::once(Domain &domain) const
{
    master->delegate->once(domain);
}
void P1D::WorkerDelegate::pass(Domain const &, PartSpecies &sp) const
{
    PartBucket L, R;
    master->delegate->partition(sp, L, R);
    //
    comm.recv<0>(master->comm.rank).unpack([&L, &R](auto payload) {
        payload.first->swap(L);
        payload.second->swap(R);
    });
    //
    sp.bucket.insert(sp.bucket.cend(), L.cbegin(), L.cend());
    sp.bucket.insert(sp.bucket.cend(), R.cbegin(), R.cend());
}
void P1D::WorkerDelegate::pass(Domain const &, ColdSpecies &sp) const
{
    recv_from_master(sp.mom0_full);
    recv_from_master(sp.mom1_full);
}
void P1D::WorkerDelegate::pass(Domain const &, BField &bfield) const
{
    recv_from_master(bfield);
}
void P1D::WorkerDelegate::pass(Domain const &, EField &efield) const
{
    recv_from_master(efield);
}
void P1D::WorkerDelegate::pass(Domain const &, Current &current) const
{
    recv_from_master(current);
}
void P1D::WorkerDelegate::gather(Domain const &, Current &current) const
{
    reduce_to_master(current);
    recv_from_master(current);
}
void P1D::WorkerDelegate::gather(Domain const &, PartSpecies &sp) const
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

template <class T, long N> void P1D::WorkerDelegate::recv_from_master(GridQ<T, N> &buffer) const
{
    comm.recv<GridQ<T, N> const *>(master->comm.rank)
        .unpack([&buffer](auto payload) noexcept(noexcept(buffer = buffer)) {
            buffer = *payload;
        });
}
template <class T, long N>
void P1D::WorkerDelegate::reduce_to_master(GridQ<T, N> const &payload) const
{
    comm.send(&payload, master->comm.rank).wait(); // must wait for delievery receipt
}
