/*
 * Copyright (c) 2022-2023, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include "Config.h"
#include "Predefined.h"
#include "UTL/Range.h"
#include <cmath>
#include <functional>
#include <iterator>
#include <map>
#include <optional>
#include <type_traits>

LIBPIC_NAMESPACE_BEGIN(1)
namespace {
template <class F>
[[nodiscard]] auto init_inverse_function_table(Range const &f_extent, Range const &x_extent, F f_of_x) -> std::map<Real, Real>
{
    static_assert(std::is_invocable_r_v<Real, F, Real>);
    std::map<Real, Real> table;
    table.insert_or_assign(end(table), f_extent.min(), x_extent.min());
    constexpr long n_samples    = 50000;
    constexpr long n_subsamples = 100;
    auto const     df           = f_extent.len / n_samples;
    auto const     dx           = x_extent.len / (n_samples * n_subsamples);
    Real           x            = x_extent.min();
    Real           f_current    = std::invoke(f_of_x, x);
    for (long i = 1; i < n_samples; ++i) {
        Real const f_target = Real(i) * df + f_extent.min();
        while (f_current < f_target)
            f_current = std::invoke(f_of_x, x += dx);
        table.insert_or_assign(end(table), f_current, x);
    }
    table.insert_or_assign(end(table), f_extent.max(), x_extent.max());
    return table;
}
[[nodiscard]] auto linear_interp(std::map<Real, Real> const &table, Real const x) noexcept -> std::optional<Real>
{
    auto const ub = table.upper_bound(x);
    if (ub == end(table) || ub == begin(table))
        return {};

    auto const &[x_min, y_min] = *std::prev(ub);
    auto const &[x_max, y_max] = *ub;
    return (y_min * (x_max - x) + y_max * (x - x_min)) / (x_max - x_min);
}

// specific to partial shell
template <bool is_small_a>
[[nodiscard]] auto int_cos_zeta(unsigned const zeta, Real const a, Real const x) noexcept -> Real
{
    constexpr auto one_sixth = 0.1666666666666666666666666666666666666667;
    auto const     ax        = a * x;
    if (zeta == 0) {
        return x;
    } else if (zeta == 1) {
        if constexpr (is_small_a) {
            return (1 - ax * ax * one_sixth) * x;
        } else {
            return std::sin(ax) / a;
        }
    } else {
        auto const addendum = int_cos_zeta<is_small_a>(zeta - 2, a, x) * Real(zeta - 1) / zeta;
        if constexpr (is_small_a) {
            return (1 + ax * ax * one_sixth * (2 - 3 * Real(zeta))) * (x / zeta) + addendum;
        } else {
            return std::pow(std::cos(ax), zeta - 1) * std::sin(ax) / (zeta * a) + addendum;
        }
    }
}
} // namespace
LIBPIC_NAMESPACE_END(1)
