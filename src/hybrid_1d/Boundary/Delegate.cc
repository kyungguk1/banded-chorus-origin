/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "Delegate.h"

#include <algorithm>
#include <random>

// MARK: Interface
//
// void H1D::Delegate::once(Domain &domain)
//{
//    std::mt19937 g{123};
//    std::uniform_real_distribution<> d{-1, 1};
//    for (Vector &v : domain.bfield) {
//        v.x += d(g) * Debug::initial_bfield_noise_amplitude;
//        v.y += d(g) * Debug::initial_bfield_noise_amplitude;
//        v.z += d(g) * Debug::initial_bfield_noise_amplitude;
//    }
//}

void H1D::Delegate::partition(PartSpecies &sp, PartBucket &L_bucket, PartBucket &R_bucket) const
{
    // note that particle position is already normalized by the grid size

    // group particles that have crossed left boundaries
    //
    auto L_it = std::partition(sp.bucket.begin(), sp.bucket.end(),
                               [LB = 0.0](Particle const &ptl) noexcept -> bool {
                                   return ptl.pos_x >= LB;
                               });
    L_bucket.insert(L_bucket.cend(), L_it, sp.bucket.end());
    sp.bucket.erase(L_it, sp.bucket.end());

    // group particles that have crossed right boundaries
    //
    auto R_it
        = std::partition(sp.bucket.begin(), sp.bucket.end(),
                         [RB = sp.params.domain_extent.len](Particle const &ptl) noexcept -> bool {
                             return ptl.pos_x < RB;
                         });
    R_bucket.insert(R_bucket.cend(), R_it, sp.bucket.end());
    sp.bucket.erase(R_it, sp.bucket.end());
}
void H1D::Delegate::pass(Domain const &domain, PartBucket &L_bucket, PartBucket &R_bucket) const
{
    // note that particle position is already normalized by the grid size

    // simulation domain size
    //
    Real const Lx = domain.params.domain_extent.len;

    for (Particle &ptl : L_bucket) { // crossed left boundary; wrap around to the rightmost cell
        ptl.pos_x += Lx;
    }
    for (Particle &ptl : R_bucket) { // crossed right boundary; wrap around to the leftmost cell
        ptl.pos_x -= Lx;
    }

    using std::swap;
    swap(L_bucket, R_bucket);
}
void H1D::Delegate::pass(Domain const &domain, PartSpecies &sp) const
{
    auto &[L, R] = buckets.cleared(); // be careful not to access it from multiple threads
                                      // be sure to clear the contents before use
    partition(sp, L, R);
    pass(domain, L, R);
    sp.bucket.insert(sp.bucket.cend(), L.cbegin(), L.cend());
    sp.bucket.insert(sp.bucket.cend(), R.cbegin(), R.cend());
}
// void H1D::Delegate::pass(Domain const&, BField &bfield) const
//{
//    if constexpr (Debug::zero_out_electromagnetic_field) {
//        bfield.fill(bfield.geomtr.B0);
//    }
//    pass(bfield);
//}
// void H1D::Delegate::pass(Domain const&, EField &efield) const
//{
//    if constexpr (Debug::zero_out_electromagnetic_field) {
//        efield.fill(Vector{});
//    }
//    pass(efield);
//}
// void H1D::Delegate::pass(Domain const&, Charge &charge) const
//{
//    pass(charge);
//}
// void H1D::Delegate::pass(Domain const&, Current &current) const
//{
//    pass(current);
//}
// void H1D::Delegate::gather(Domain const&, Charge &charge) const
//{
//    gather(charge);
//}
// void H1D::Delegate::gather(Domain const&, Current &current) const
//{
//    gather(current);
//}
// void H1D::Delegate::gather(Domain const&, PartSpecies &sp) const
//{
//    gather(sp.moment<0>());
//    gather(sp.moment<1>());
//    gather(sp.moment<2>());
//}

// MARK: Implementation
//
template <class T, long N> void H1D::Delegate::pass(GridQ<T, N> &A)
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
template <class T, long N> void H1D::Delegate::gather(GridQ<T, N> &A)
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
