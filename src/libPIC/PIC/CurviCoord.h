/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC/Config.h>

#include <ostream>
#include <type_traits>

LIBPIC_BEGIN_NAMESPACE
struct CurviCoord {
    using Real = double;

    Real q1{};

    // constructors
    //
    constexpr CurviCoord() noexcept = default;
    constexpr explicit CurviCoord(Real q1) noexcept
    : q1{ q1 } {}

    // compound operations
    //
    constexpr CurviCoord &operator+=(CurviCoord const &o) noexcept
    {
        q1 += o.q1;
        return *this;
    }
    constexpr CurviCoord &operator-=(CurviCoord const &o) noexcept
    {
        q1 -= o.q1;
        return *this;
    }
    constexpr CurviCoord &operator*=(CurviCoord const &o) noexcept
    {
        q1 *= o.q1;
        return *this;
    }
    constexpr CurviCoord &operator/=(CurviCoord const &o) noexcept
    {
        q1 /= o.q1;
        return *this;
    }

    // unary operations
    //
    [[nodiscard]] friend constexpr decltype(auto) operator+(CurviCoord const &s) noexcept
    {
        return s;
    }
    [[nodiscard]] friend constexpr auto operator-(CurviCoord const &s) noexcept
    {
        return CurviCoord{ -s.q1 };
    }

    // binary operations
    //
    [[nodiscard]] friend constexpr auto operator+(CurviCoord a, CurviCoord const &b) noexcept
    {
        a += b;
        return a;
    }
    [[nodiscard]] friend constexpr auto operator-(CurviCoord a, CurviCoord const &b) noexcept
    {
        a -= b;
        return a;
    }
    [[nodiscard]] friend constexpr auto operator*(CurviCoord a, CurviCoord const &b) noexcept
    {
        a *= b;
        return a;
    }
    [[nodiscard]] friend constexpr auto operator/(CurviCoord a, CurviCoord const &b) noexcept
    {
        a /= b;
        return a;
    }

    // pretty print
    //
    template <class CharT, class Traits>
    friend decltype(auto) operator<<(std::basic_ostream<CharT, Traits> &os, CurviCoord const &s)
    {
        return os << '{' << s.q1 << '}';
    }
};

static_assert(std::is_standard_layout_v<CurviCoord>);
LIBPIC_END_NAMESPACE
