/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "Delegate.h"

#include <algorithm>
#include <cmath>
#include <random>

PIC1D_BEGIN_NAMESPACE
void Delegate::partition(PartSpecies &sp, PartBucket &L_bucket, PartBucket &R_bucket) const
{
    // group particles that have crossed left boundary
    auto L_it
        = std::partition(sp.bucket.begin(), sp.bucket.end(),
                         [LB = sp.params.full_grid_subdomain_extent.min()](Particle const &ptl) noexcept -> bool {
                             return ptl.pos.q1 >= LB;
                         });
    L_bucket.insert(L_bucket.cend(), L_it, sp.bucket.end());
    sp.bucket.erase(L_it, sp.bucket.end());

    // group particles that have crossed right boundary
    auto R_it
        = std::partition(sp.bucket.begin(), sp.bucket.end(),
                         [RB = sp.params.full_grid_subdomain_extent.max()](Particle const &ptl) noexcept -> bool {
                             return ptl.pos.q1 < RB;
                         });
    R_bucket.insert(R_bucket.cend(), R_it, sp.bucket.end());
    sp.bucket.erase(R_it, sp.bucket.end());
}
void Delegate::boundary_pass(Domain const &domain, PartBucket &L_bucket, PartBucket &R_bucket) const
{
    switch (domain.params.particle_boundary_condition) {
        case BC::periodic:
            periodic_particle_pass(domain.params, L_bucket, R_bucket);
            break;
        case BC::reflecting:
            reflecting_particle_pass(domain.params, L_bucket, R_bucket);
            break;
    }

    using std::swap;
    swap(L_bucket, R_bucket);
}
void Delegate::boundary_pass(Domain const &domain, PartSpecies &sp) const
{
    auto &[L, R] = buckets.cleared(); // be careful not to access it from multiple threads
                                      // be sure to clear the contents before use
    partition(sp, L, R);
    boundary_pass(domain, L, R);
    sp.bucket.insert(sp.bucket.cend(), L.cbegin(), L.cend());
    sp.bucket.insert(sp.bucket.cend(), R.cbegin(), R.cend());
}

template <class T, long N>
void Delegate::pass(Grid<T, N, Pad> &A)
{
    // fill ghost cells
    //
    for (long p = 0, m = -1; p < Pad; ++p, --m) {
        // put left boundary value to right ghost cell
        A.end()[p] = A[p];
        // put right boundary value to left ghost cell
        A[m] = A.end()[m];
    }
}
template <class T, long N>
void Delegate::gather(Grid<T, N, Pad> &A)
{
    // gather moments at ghost cells
    //
    for (long p = Pad - 1, m = -Pad; m < 0; --p, ++m) {
        // add right ghost cell value to left boundary
        A[p] += A.end()[p];
        // add left ghost cell value to right boundary
        A.end()[m] += A[m];
    }
}

void Delegate::periodic_particle_pass(ParamSet const &params, PartBucket &L_bucket, PartBucket &R_bucket)
{
    // adjust coordinates
    auto const extent = params.full_grid_whole_domain_extent;
    for (Particle &ptl : L_bucket) {
        // crossed left boundary; wrap around to the rightmost cell
        if (!extent.is_member(ptl.pos.q1))
            ptl.pos.q1 += extent.len;
    }
    for (Particle &ptl : R_bucket) {
        // crossed right boundary; wrap around to the leftmost cell
        if (!extent.is_member(ptl.pos.q1))
            ptl.pos.q1 -= extent.len;
    }
}
void Delegate::reflecting_particle_pass(ParamSet const &params, PartBucket &L_bucket, PartBucket &R_bucket)
{
    // adjust coordinates and flip velocity vector
    auto const extent = params.full_grid_whole_domain_extent;
    for (Particle &ptl : L_bucket) {
        if (!extent.is_member(ptl.pos.q1)) {
            // crossed left boundary; put back into the leftmost cell
            ptl.pos.q1 += 2 * (extent.min() - ptl.pos.q1);
            // velocity flip
            ptl.g_vel.x *= -1;
            // gyro-phase randomization
            if constexpr (ParamSet::should_randomize_gyrophase_of_reflecting_particles)
                ptl.g_vel = randomize_gyrophase(ptl.g_vel);
        }
    }
    for (Particle &ptl : R_bucket) {
        if (!extent.is_member(ptl.pos.q1)) {
            // crossed right boundary; put back into the rightmost cell
            ptl.pos.q1 += 2 * (extent.max() - ptl.pos.q1);
            // velocity flip
            ptl.g_vel.x *= -1;
            // gyro-phase randomization
            if constexpr (ParamSet::should_randomize_gyrophase_of_reflecting_particles)
                ptl.g_vel = randomize_gyrophase(ptl.g_vel);
        }
    }
}
auto Delegate::randomize_gyrophase(Vector const &v) -> Vector
{
    thread_local static std::mt19937        rng{ std::random_device{}() };
    static std::uniform_real_distribution<> dist{ -M_PI, M_PI };

    auto const phase = dist(rng);
    auto const cos   = std::cos(phase);
    auto const sin   = std::sin(phase);

    return { v.x, v.y * cos - v.z * sin, v.y * sin + v.z * cos };
}
PIC1D_END_NAMESPACE
