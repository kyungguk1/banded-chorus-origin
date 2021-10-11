/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "Delegate.h"

#include <algorithm>

HYBRID1D_BEGIN_NAMESPACE
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
    // simulation domain extent
    auto const extent = domain.params.full_grid_whole_domain_extent;

    // periodic boundary
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
HYBRID1D_END_NAMESPACE
