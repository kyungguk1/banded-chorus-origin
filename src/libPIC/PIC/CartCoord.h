/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC/Config.h>
#include <PIC/Scalar.h>

#include <ostream>
#include <type_traits>

LIBPIC_BEGIN_NAMESPACE
class CartCoord : private Scalar {
public:
    // value access
    //
    using Scalar::operator Real;
    using Scalar::operator*;
    using Scalar::operator();

    // constructors
    //
    using Scalar::Scalar;

    // compound operations
    //
    constexpr decltype(auto) operator+=(CartCoord const &o) noexcept
    {
        return static_cast<CartCoord &>(Scalar::operator+=(o));
    }
    constexpr decltype(auto) operator-=(CartCoord const &o) noexcept
    {
        return static_cast<CartCoord &>(Scalar::operator-=(o));
    }
    constexpr decltype(auto) operator*=(CartCoord const &o) noexcept
    {
        return static_cast<CartCoord &>(Scalar::operator*=(o));
    }
    constexpr decltype(auto) operator/=(CartCoord const &o) noexcept
    {
        return static_cast<CartCoord &>(Scalar::operator/=(o));
    }

    // unary operations
    //
    [[nodiscard]] friend constexpr CartCoord const &operator+(CartCoord const &s) noexcept
    {
        return s;
    }
    [[nodiscard]] friend constexpr CartCoord operator-(CartCoord const &s) noexcept
    {
        return -Real{ s };
    }

    // binary operations
    //
    [[nodiscard]] friend constexpr CartCoord operator+(CartCoord a, CartCoord const &b) noexcept
    {
        return a += b;
    }
    [[nodiscard]] friend constexpr CartCoord operator-(CartCoord a, CartCoord const &b) noexcept
    {
        return a -= b;
    }
    [[nodiscard]] friend constexpr CartCoord operator*(CartCoord a, CartCoord const &b) noexcept
    {
        return a *= b;
    }
    [[nodiscard]] friend constexpr CartCoord operator/(CartCoord a, CartCoord const &b) noexcept
    {
        return a /= b;
    }

    // pretty print
    //
    template <class CharT, class Traits>
    friend decltype(auto) operator<<(std::basic_ostream<CharT, Traits> &os, CartCoord const &s)
    {
        return os << static_cast<Scalar const &>(s);
    }
};

static_assert(std::is_standard_layout_v<CartCoord>);
LIBPIC_END_NAMESPACE
