/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include "../Boundary/Delegate.h"
#include "Domain.h"

HYBRID1D_BEGIN_NAMESPACE
namespace {
template <class T, long N>
auto &operator+=(Grid<T, N, Pad> &lhs, Grid<T, N, Pad> const &rhs) noexcept
{
    auto rhs_first = rhs.dead_begin(), rhs_last = rhs.dead_end();
    auto lhs_first = lhs.dead_begin();
    while (rhs_first != rhs_last) {
        *lhs_first++ += *rhs_first++;
    }
    return lhs;
}
//
template <class T, long N> auto &operator*=(Grid<T, N, Pad> &lhs, T const rhs) noexcept
{
    auto first = lhs.dead_begin(), last = lhs.dead_end();
    while (first != last) {
        *first++ *= rhs;
    }
    return lhs;
}
} // namespace

template <class Species>
auto Domain::collect_smooth(Charge &rho, Species const &sp) const -> Charge const &
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
        delegate->pass(*this, rho);
        rho.smooth();
    }
    //
    delegate->pass(*this, rho);
    return rho;
}
template <class Species>
auto Domain::collect_smooth(Current &J, Species const &sp) const -> Current const &
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
        delegate->pass(*this, J);
        J.smooth();
    }
    //
    delegate->pass(*this, J);
    return J;
}
HYBRID1D_END_NAMESPACE
