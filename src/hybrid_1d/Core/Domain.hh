/*
 * Copyright (c) 2019, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#ifndef Domain_hh
#define Domain_hh

#include "../Boundary/Delegate.h"
#include "Domain.h"

template <class Species>
auto H1D::Domain::collect_smooth(Charge &rho, Species const &sp) const -> Charge const &
{
    rho.reset();
    //
    // collect & gather rho
    //
    rho += sp;
    delegate->gather(*this, rho);
    //
    // optional smoothing
    //
    for (long i = 0; i < sp->number_of_source_smoothings; ++i) {
        delegate->pass(*this, rho), rho.smooth();
    }
    //
    return delegate->pass(*this, rho), rho;
}
template <class Species>
auto H1D::Domain::collect_smooth(Current &J, Species const &sp) const -> Current const &
{
    J.reset();
    //
    // collect & gather J
    //
    J += sp;
    delegate->gather(*this, J);
    //
    // optional smoothing
    //
    for (long i = 0; i < sp->number_of_source_smoothings; ++i) {
        delegate->pass(*this, J), J.smooth();
    }
    //
    return delegate->pass(*this, J), J;
}

namespace {
template <class T, long N>
auto &operator+=(H1D::GridQ<T, N> &lhs, H1D::GridQ<T, N> const &rhs) noexcept
{
    auto rhs_first = rhs.dead_begin(), rhs_last = rhs.dead_end();
    auto lhs_first = lhs.dead_begin();
    while (rhs_first != rhs_last) {
        *lhs_first++ += *rhs_first++;
    }
    return lhs;
}
//
template <class T, long N> auto &operator*=(H1D::GridQ<T, N> &lhs, T const rhs) noexcept
{
    auto first = lhs.dead_begin(), last = lhs.dead_end();
    while (first != last) {
        *first++ *= rhs;
    }
    return lhs;
}
} // namespace

#endif /* Domain_hh */
