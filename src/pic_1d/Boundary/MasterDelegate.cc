/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "MasterDelegate.h"

#include <algorithm>
#include <functional>
#include <iterator>
#include <utility>

PIC1D_BEGIN_NAMESPACE
namespace {
template <class T, long N>
decltype(auto) operator/=(Grid<T, N, Pad> &G, T const w) noexcept
{
    // include padding
    std::for_each(G.dead_begin(), G.dead_end(), [w](T &value_ref) noexcept {
        value_ref /= w;
    });
    return G;
}
} // namespace

MasterDelegate::~MasterDelegate()
{
}
MasterDelegate::MasterDelegate(Delegate *const delegate)
: delegate{ delegate }, all_but_master{}
{
    comm = dispatch.comm(static_cast<unsigned>(workers.size()));
    for (unsigned i = 0; i < workers.size(); ++i) {
        workers[i].master = this;
        workers[i].comm   = dispatch.comm(i);
        all_but_master.emplace_back(i);
    }
}

void MasterDelegate::setup(Domain &domain) const
{
    // distribute particles to workers
    //
    for (PartSpecies &sp : domain.part_species) {
        sp.equilibrium_macro_weight(Badge<MasterDelegate>{}) /= ParamSet::number_of_particle_parallelism;
        distribute(domain, sp);
    }

    // distribute cold species moments to workers
    //
    for (ColdSpecies &sp : domain.cold_species) {
        // evenly split moments
        auto const divisor = Real(workers.size() + 1);
        sp.mom0_full /= Scalar{ divisor };
        sp.mom1_full /= Vector{ divisor };

        // pass along
        distribute(domain, sp);
    }
}
void MasterDelegate::distribute(Domain const &, PartSpecies &sp) const
{
    // distribute particles to workers
    //
    std::vector<decltype(sp.bucket)> payloads;
    payloads.reserve(all_but_master.size());
    auto const chunk = static_cast<long>(sp.bucket.size() / (workers.size() + 1));
    for ([[maybe_unused]] rank_t const &rank : all_but_master) { // master excluded
        auto const last  = end(sp.bucket);
        auto const first = std::prev(last, chunk);
        payloads.emplace_back(std::make_move_iterator(first), std::make_move_iterator(last));
        sp.bucket.erase(first, last);
    }
    auto tks = comm.scatter(std::move(payloads), all_but_master);
    std::for_each(std::make_move_iterator(begin(tks)), std::make_move_iterator(end(tks)),
                  std::mem_fn(&ticket_t::wait));
}
void MasterDelegate::distribute(Domain const &, ColdSpecies &sp) const
{
    // distribute cold species moments to workers
    //
    broadcast_to_workers(sp.mom0_full);
    broadcast_to_workers(sp.mom1_full);
}

void MasterDelegate::teardown(Domain &domain) const
{
    // collect particles from workers
    //
    for (PartSpecies &sp : domain.part_species) {
        collect(domain, sp);
        sp.equilibrium_macro_weight(Badge<MasterDelegate>{}) *= ParamSet::number_of_particle_parallelism;
    }

    // collect cold species from workers
    //
    for (ColdSpecies &sp : domain.cold_species) {
        // moments are automatically accumulated
        collect(domain, sp);
    }
}
void MasterDelegate::collect(Domain const &, PartSpecies &sp) const
{
    // gather particles from workers
    //
    using PartBucket = decltype(sp.bucket);
    comm.for_each<PartBucket>(
        all_but_master,
        [](PartBucket payload, PartBucket &bucket) {
            std::move(begin(payload), end(payload), std::back_inserter(bucket));
        },
        sp.bucket);
}
void MasterDelegate::collect(Domain const &, ColdSpecies &sp) const
{
    collect_from_workers(sp.mom0_full);
    collect_from_workers(sp.mom1_full);
}

void MasterDelegate::prologue(Domain const &domain, long const i) const
{
    delegate->prologue(domain, i);
}
void MasterDelegate::epilogue(Domain const &domain, long const i) const
{
    delegate->epilogue(domain, i);
}
void MasterDelegate::once(Domain &domain) const
{
    delegate->once(domain);
}
void MasterDelegate::boundary_pass(Domain const &domain, PartSpecies &sp) const
{
    auto &[L, R] = buckets.cleared(); // be careful not to access it from multiple threads
                                      // be sure to clear the contents before use
    delegate->partition(sp, L, R);
    //
    delegate->boundary_pass(domain, L, R);
    for (auto const &worker : workers) {
        comm.send(std::make_pair(&L, &R), worker.comm.rank).wait();
        delegate->boundary_pass(domain, L, R);
    }
    //
    sp.bucket.insert(sp.bucket.cend(), L.cbegin(), L.cend());
    sp.bucket.insert(sp.bucket.cend(), R.cbegin(), R.cend());
}
void MasterDelegate::boundary_pass(Domain const &domain, ColdSpecies &sp) const
{
    delegate->boundary_pass(domain, sp);
    broadcast_to_workers(sp.mom0_full);
    broadcast_to_workers(sp.mom1_full);
}
void MasterDelegate::boundary_pass(Domain const &domain, BField &bfield) const
{
    delegate->boundary_pass(domain, bfield);
    broadcast_to_workers(bfield);
}
void MasterDelegate::boundary_pass(Domain const &domain, EField &efield) const
{
    delegate->boundary_pass(domain, efield);
    broadcast_to_workers(efield);
}
void MasterDelegate::boundary_pass(Domain const &domain, Current &current) const
{
    delegate->boundary_pass(domain, current);
    broadcast_to_workers(current);
}
void MasterDelegate::boundary_gather(Domain const &domain, Current &current) const
{
    collect_from_workers(current);
    delegate->boundary_gather(domain, current);
    broadcast_to_workers(current);
}
void MasterDelegate::boundary_gather(Domain const &domain, Species &sp) const
{
    {
        collect_from_workers(sp.moment<0>());
        collect_from_workers(sp.moment<1>());
        collect_from_workers(sp.moment<2>());
    }
    delegate->boundary_gather(domain, sp);
    {
        broadcast_to_workers(sp.moment<0>());
        broadcast_to_workers(sp.moment<1>());
        broadcast_to_workers(sp.moment<2>());
    }
}

namespace {
template <class T, long N, class U>
decltype(auto) operator/=(Grid<T, N, Pad> &lhs, U const w) noexcept
{ // include padding
    for (auto it = lhs.dead_begin(), end = lhs.dead_end(); it != end; ++it) {
        *it /= w;
    }
    return lhs;
}
template <class T, long N>
decltype(auto) operator+=(Grid<T, N, Pad> &lhs, Grid<T, N, Pad> const &rhs) noexcept
{
    auto lhs_first = lhs.dead_begin(), lhs_last = lhs.dead_end();
    auto rhs_first = rhs.dead_begin();
    while (lhs_first != lhs_last) {
        *lhs_first++ += *rhs_first++;
    }
    return lhs;
}
} // namespace
template <class T, long N>
void MasterDelegate::broadcast_to_workers(Grid<T, N, Pad> const &payload) const
{
    auto tks = comm.bcast(&payload, all_but_master);
    std::for_each(std::make_move_iterator(begin(tks)), std::make_move_iterator(end(tks)),
                  std::mem_fn(&ticket_t::wait));
}
template <class T, long N>
void MasterDelegate::collect_from_workers(Grid<T, N, Pad> &buffer) const
{
    // the first worker will collect all workers'
    //
    comm.for_each<Grid<T, N, Pad> const *>(
        all_but_master,
        [](auto payload, Grid<T, N, Pad> &buffer) {
            buffer += *payload;
        },
        buffer);
}
PIC1D_END_NAMESPACE
