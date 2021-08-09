/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "SubdomainDelegate.h"

#include <algorithm>
#include <iterator>
#include <random>
#include <stdexcept>
#include <utility>

// MARK:- P1D::SubdomainDelegate
//
P1D::SubdomainDelegate::SubdomainDelegate(parallel::mpi::Comm _comm) : comm{ std::move(_comm) }
{
    if (!comm->operator bool())
        throw std::invalid_argument{ __PRETTY_FUNCTION__ };

    int const size = comm.size();
    int const rank = comm->rank();
    left_          = rank_t{ (size + rank - 1) % size };
    right          = rank_t{ (size + rank + 1) % size };
}

void P1D::SubdomainDelegate::once(Domain &domain) const
{
    std::mt19937                     g{ 494983U + static_cast<unsigned>(comm->rank()) };
    std::uniform_real_distribution<> d{ -1, 1 };
    for (Vector &v : domain.efield) {
        v.x += d(g) * Debug::initial_efield_noise_amplitude;
        v.y += d(g) * Debug::initial_efield_noise_amplitude;
        v.z += d(g) * Debug::initial_efield_noise_amplitude;
    }
}

void P1D::SubdomainDelegate::pass(Domain const &, ColdSpecies &sp) const
{
    pass(sp.mom0_full);
    pass(sp.mom1_full);
}
void P1D::SubdomainDelegate::pass(Domain const &, BField &bfield) const
{
    if constexpr (Debug::zero_out_electromagnetic_field) {
        bfield.fill(bfield.geomtr.B0);
    } else if constexpr (Input::is_electrostatic) { // zero-out transverse components
        for (Vector &v : bfield) {
            v.y = bfield.geomtr.B0.y;
            v.z = bfield.geomtr.B0.z;
        }
    }
    pass(bfield);
}
void P1D::SubdomainDelegate::pass(Domain const &, EField &efield) const
{
    if constexpr (Debug::zero_out_electromagnetic_field) {
        efield.fill(Vector{});
    } else if constexpr (Input::is_electrostatic) { // zero-out transverse components
        for (Vector &v : efield) {
            v.y = v.z = 0;
        }
    }
    pass(efield);
}
void P1D::SubdomainDelegate::pass(Domain const &, Current &current) const
{
    pass(current);
}
void P1D::SubdomainDelegate::gather(Domain const &, Current &current) const
{
    gather(current);
}
void P1D::SubdomainDelegate::gather(Domain const &, PartSpecies &sp) const
{
    gather(sp.moment<0>());
    gather(sp.moment<1>());
    gather(sp.moment<2>());
}
void P1D::SubdomainDelegate::pass(Domain const &domain, PartBucket &L_bucket,
                                  PartBucket &R_bucket) const
{
    // pass across boundaries
    // send-recv pair order is important
    // e.g., if send-left is first, recv-right should appear first.
    //
    constexpr parallel::mpi::Tag tag1{ 1 }, tag2{ 2 };
    {
        auto tk1 = comm.ibsend(std::move(L_bucket), { left_, tag1 });
        auto tk2 = comm.ibsend(std::move(R_bucket), { right, tag2 });
        {
            L_bucket = comm.recv<3>({}, { right, tag1 });
            R_bucket = comm.recv<3>({}, { left_, tag2 });
        }
        std::move(tk1).wait();
        std::move(tk2).wait();
    }

    // adjust coordinates
    //
    Delegate::pass(domain, L_bucket, R_bucket);
}
template <class T, long Mx> void P1D::SubdomainDelegate::pass(GridQ<T, Mx> &grid) const
{
    // pass across boundaries
    // send-recv pair order is important
    // e.g., if send-left is first, recv-right should appear first.

    constexpr parallel::mpi::Tag tag1{ 1 }, tag2{ 2 };
    if constexpr (Mx >= Pad) {
        auto tk_left_ = comm.issend<T>(grid.begin(), std::next(grid.begin(), Pad), { left_, tag1 });
        auto tk_right = comm.issend<T>(std::prev(grid.end(), Pad), grid.end(), { right, tag2 });
        {
            comm.recv<T>(grid.end(), std::next(grid.end(), Pad), { right, tag1 });
            comm.recv<T>(std::prev(grid.begin(), Pad), grid.begin(), { left_, tag2 });
        }
        std::move(tk_left_).wait();
        std::move(tk_right).wait();
    } else {
        // from inside out
        //
        for (long b = 0, e = -1; b < Pad; ++b, --e) {
            auto tk_left_ = comm.issend<T>(grid.begin()[b], { left_, tag1 });
            auto tk_right = comm.issend<T>(grid.end()[e], { right, tag2 });
            {
                grid.end()[b]   = comm.recv<T>({ right, tag1 });
                grid.begin()[e] = comm.recv<T>({ left_, tag2 });
            }
            std::move(tk_left_).wait();
            std::move(tk_right).wait();
        }
    }
}
template <class T, long Mx> void P1D::SubdomainDelegate::gather(GridQ<T, Mx> &grid) const
{
    // pass across boundaries
    // send-recv pair order is important
    // e.g., if send-left is first, recv-right should appear first.

    constexpr parallel::mpi::Tag tag1{ 1 }, tag2{ 2 };
    if constexpr (Mx >= Pad) {
        auto accum = [](auto payload, auto *first, auto *last) {
            std::transform(first, last, begin(payload), first, std::plus{});
        };

        auto tk_left_ = comm.issend<T>(std::prev(grid.begin(), Pad), grid.begin(), { left_, tag1 });
        auto tk_right = comm.issend<T>(grid.end(), std::next(grid.end(), Pad), { right, tag2 });
        {
            comm.recv<T>({}, { right, tag1 }).unpack(accum, std::prev(grid.end(), Pad), grid.end());
            comm.recv<T>({}, { left_, tag2 })
                .unpack(accum, grid.begin(), std::next(grid.begin(), Pad));
        }
        std::move(tk_left_).wait();
        std::move(tk_right).wait();
    } else {
        // from outside in
        //
        for (long b = -Pad, e = Pad - 1; b < 0; ++b, --e) {
            auto tk_left_ = comm.issend<T>(grid.begin()[b], { left_, tag1 });
            auto tk_right = comm.issend<T>(grid.end()[e], { right, tag2 });
            {
                grid.end()[b] += comm.recv<T>({ right, tag1 });
                grid.begin()[e] += comm.recv<T>({ left_, tag2 });
            }
            std::move(tk_left_).wait();
            std::move(tk_right).wait();
        }
    }
}
