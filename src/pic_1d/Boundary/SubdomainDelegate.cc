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

PIC1D_BEGIN_NAMESPACE
SubdomainDelegate::SubdomainDelegate(parallel::mpi::Comm _comm)
: comm{ std::move(_comm) }
{
    if (!comm->operator bool())
        throw std::invalid_argument{ __PRETTY_FUNCTION__ };

    int const size = comm.size();
    int const rank = comm->rank();
    left_          = rank_t{ (size + rank - 1) % size };
    right          = rank_t{ (size + rank + 1) % size };
}

void SubdomainDelegate::once(Domain &domain) const
{
    std::mt19937                     g{ 494983U + static_cast<unsigned>(comm->rank()) };
    std::uniform_real_distribution<> d{ -1, 1 };
    for (Vector &v : domain.efield) {
        v.x += d(g) * Debug::initial_efield_noise_amplitude;
        v.y += d(g) * Debug::initial_efield_noise_amplitude;
        v.z += d(g) * Debug::initial_efield_noise_amplitude;
    }
}

void SubdomainDelegate::boundary_pass(Domain const &, ColdSpecies &sp) const
{
    mpi_pass(sp.mom0_full);
    mpi_pass(sp.mom1_full);
}
void SubdomainDelegate::boundary_pass(Domain const &, BField &bfield) const
{
    if constexpr (Debug::zero_out_electromagnetic_field) {
        bfield.fill(Vector{});
    } else if constexpr (Input::is_electrostatic) { // zero-out transverse components
        for (Vector &v : bfield) {
            v.y = 0;
            v.z = 0;
        }
    }
    mask(bfield.params, bfield);
    mpi_pass(bfield);
}
void SubdomainDelegate::boundary_pass(Domain const &, EField &efield) const
{
    if constexpr (Debug::zero_out_electromagnetic_field) {
        efield.fill(Vector{});
    } else if constexpr (Input::is_electrostatic) { // zero-out transverse components
        for (Vector &v : efield) {
            v.y = 0;
            v.z = 0;
        }
    }
    mask(efield.params, efield);
    mpi_pass(efield);
}
template <class T, long Mx>
void SubdomainDelegate::mask(ParamSet const &params, Grid<T, Mx, Pad> &grid) const
{
    auto const left_offset  = params.full_grid_subdomain_extent.min() - params.full_grid_whole_domain_extent.min();
    auto const right_offset = params.full_grid_whole_domain_extent.max() - params.full_grid_subdomain_extent.max();
    for (long i = 0, first = 0, last = Mx - 1; i < Mx; ++i) {
        grid[first++] *= params.masking_function(left_offset + i);
        grid[last--] *= params.masking_function(right_offset + i);
    }
}
template <class T, long Mx>
void SubdomainDelegate::mpi_pass(Grid<T, Mx, Pad> &grid) const
{
    // pass across boundaries
    // send-recv pair order is important
    // e.g., if send-left is first, recv-right should appear first.

    constexpr parallel::mpi::Tag tag1{ 1 };
    constexpr parallel::mpi::Tag tag2{ 2 };
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

void SubdomainDelegate::boundary_pass(Domain const &, Current &current) const
{
    mpi_pass(current);
}
void SubdomainDelegate::boundary_gather(Domain const &, Current &current) const
{
    moment_gather(current.params, current);

    mpi_gather(current);
}
void SubdomainDelegate::boundary_gather(Domain const &, Species &sp) const
{
    moment_gather(sp.params, sp.moment<0>());
    moment_gather(sp.params, sp.moment<1>());
    moment_gather(sp.params, sp.moment<2>());

    mpi_gather(sp.moment<0>());
    mpi_gather(sp.moment<1>());
    mpi_gather(sp.moment<2>());
}
template <class T, long Mx>
void SubdomainDelegate::moment_gather(ParamSet const &params, Grid<T, Mx, Pad> &grid) const
{
    switch (params.particle_boundary_condition) {
        case BC::periodic:
            // do nothing
            break;
        case BC::reflecting: {
            if (is_leftmost_subdomain()) {
                // at the leftmost subdomain, moments are accumulated at the first grid point
                for (long i = -Pad; i < 0; ++i) {
                    grid[i + 1] += std::exchange(grid[i], T{});
                }
            }
            if (is_rightmost_subdomain()) {
                // at the rightmost subdomain, moments are accumulated at the last grid point
                for (long i = Mx + Pad - 1; i >= Mx; --i) {
                    grid[i - 1] += std::exchange(grid[i], T{});
                }
            }
            break;
        }
    }
}
template <class T, long Mx>
void SubdomainDelegate::mpi_gather(Grid<T, Mx, Pad> &grid) const
{
    // pass across boundaries
    // send-recv pair order is important
    // e.g., if send-left is first, recv-right should appear first.

    constexpr parallel::mpi::Tag tag1{ 1 };
    constexpr parallel::mpi::Tag tag2{ 2 };
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

void SubdomainDelegate::boundary_pass(Domain const &domain, PartBucket &L_bucket, PartBucket &R_bucket) const
{
    // pass particles across subdomain boundaries
    //
    switch (domain.params.particle_boundary_condition) {
        case BC::periodic:
            periodic_particle_pass(domain, L_bucket, R_bucket);
            break;
        case BC::reflecting:
            reflecting_particle_pass(domain, L_bucket, R_bucket);
            break;
    }

    // adjust coordinates and (if necessary) flip velocity vector
    //
    Delegate::boundary_pass(domain, L_bucket, R_bucket);
}
void SubdomainDelegate::periodic_particle_pass(Domain const &, PartBucket &L_bucket, PartBucket &R_bucket) const
{
    // pass particles across boundaries
    //
    mpi_pass(L_bucket, R_bucket);
}
void SubdomainDelegate::reflecting_particle_pass(Domain const &, PartBucket &L_bucket, PartBucket &R_bucket) const
{
    // hijack boundary-crossing particles
    //
    PartBucket L_hijacked;
    if (is_leftmost_subdomain()) {
        L_hijacked = std::move(L_bucket);
        L_bucket.clear();
    } else {
        L_hijacked.clear();
    }

    PartBucket R_hijacked{};
    if (is_rightmost_subdomain()) {
        R_hijacked = std::move(R_bucket);
        R_bucket.clear();
    } else {
        R_hijacked.clear();
    }

    // pass remaining particles across subdomain boundaries
    //
    mpi_pass(L_bucket, R_bucket);

    // put back the hijacked particles
    //
    L_bucket.insert(end(L_bucket), std::make_move_iterator(begin(L_hijacked)), std::make_move_iterator(end(L_hijacked)));
    R_bucket.insert(end(R_bucket), std::make_move_iterator(begin(R_hijacked)), std::make_move_iterator(end(R_hijacked)));
}
void SubdomainDelegate::mpi_pass(PartBucket &L_bucket, PartBucket &R_bucket) const
{
    // pass particles across boundaries
    // send-recv pair order is important
    // e.g., if send-left is first, recv-right should appear first.
    //
    constexpr parallel::mpi::Tag tag1{ 1 };
    constexpr parallel::mpi::Tag tag2{ 2 };
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
}
PIC1D_END_NAMESPACE
