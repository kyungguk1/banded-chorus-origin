/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include "PIC/BitReversedPattern.h"
#include "PIC/xoroshiro128.h"
#include <PIC/Config.h>
#include <PIC/Predefined.h>

#include <random>

LIBPIC_BEGIN_NAMESPACE
// given a rng object, returns a real number, (0, 1), following a uniform distribution
//
template <class URBG> [[nodiscard]] static Real uniform_real(URBG &g) noexcept
{
    using uniform_t                   = std::uniform_real_distribution<>;
    static constexpr Real         eps = 1e-15;
    thread_local static uniform_t uniform{ eps, 1 - eps };
    return uniform(g);
}

template <unsigned seed> [[nodiscard]] static Real uniform_mt19937() noexcept
{
    thread_local static std::mt19937 g{ seed };
    return uniform_real(g);
}
template <unsigned seed> [[nodiscard]] static Real uniform_xoroshiro128() noexcept
{
    thread_local static xoroshiro128<seed> g{};
    return uniform_real(g);
}

/// Returns a real number (0, 1) following a uniform distribution
/// \tparam seed A seed for a random number generator.
///
template <unsigned seed> [[nodiscard]] static Real uniform_real() noexcept
{ // seed must be passed as a template parameter
    return uniform_mt19937<seed>();
}

/// Returns a real number (0, 1) following a uniform distribution
/// \tparam base Base prime number for BitReversedPattern.
///
template <unsigned base> [[nodiscard]] static Real bit_reversed() noexcept
{
    static_assert(base > 0, "base has to be a positive number");
    thread_local static BitReversedPattern<base> g{};
    return uniform_real(g);
}
LIBPIC_END_NAMESPACE
