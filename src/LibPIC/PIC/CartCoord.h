/*
 * Copyright (c) 2021-2022, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC/Config.h>

#include <ostream>
#include <type_traits>

LIBPIC_NAMESPACE_BEGIN(1)
struct CartCoord {
    using Real = double;

    Real x{};

    // constructors
    //
    constexpr CartCoord() noexcept = default;
    constexpr explicit CartCoord(Real x) noexcept
    : x{ x } {}

    // compound operations
    //
    constexpr CartCoord &operator+=(CartCoord const &o) noexcept
    {
        x += o.x;
        return *this;
    }
    constexpr CartCoord &operator-=(CartCoord const &o) noexcept
    {
        x -= o.x;
        return *this;
    }
    constexpr CartCoord &operator*=(CartCoord const &o) noexcept
    {
        x *= o.x;
        return *this;
    }
    constexpr CartCoord &operator/=(CartCoord const &o) noexcept
    {
        x /= o.x;
        return *this;
    }

    // unary operations
    //
    [[nodiscard]] friend constexpr decltype(auto) operator+(CartCoord const &s) noexcept
    {
        return s;
    }
    [[nodiscard]] friend constexpr auto operator-(CartCoord const &s) noexcept
    {
        return CartCoord{ -s.x };
    }

    // binary operations
    //
    [[nodiscard]] friend constexpr auto operator+(CartCoord a, CartCoord const &b) noexcept
    {
        a += b;
        return a;
    }
    [[nodiscard]] friend constexpr auto operator-(CartCoord a, CartCoord const &b) noexcept
    {
        a -= b;
        return a;
    }
    [[nodiscard]] friend constexpr auto operator*(CartCoord a, CartCoord const &b) noexcept
    {
        a *= b;
        return a;
    }
    [[nodiscard]] friend constexpr auto operator/(CartCoord a, CartCoord const &b) noexcept
    {
        a /= b;
        return a;
    }

    // pretty print
    //
    template <class CharT, class Traits>
    friend decltype(auto) operator<<(std::basic_ostream<CharT, Traits> &os, CartCoord const &s)
    {
        return os << '{' << s.x << '}';
    }
};

static_assert(std::is_standard_layout_v<CartCoord>);
LIBPIC_NAMESPACE_END(1)
