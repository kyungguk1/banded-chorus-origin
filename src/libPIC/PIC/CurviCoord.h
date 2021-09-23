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
class CurviCoord : private Scalar {
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
    constexpr decltype(auto) operator+=(CurviCoord const &o) noexcept
    {
        return static_cast<CurviCoord &>(Scalar::operator+=(o));
    }
    constexpr decltype(auto) operator-=(CurviCoord const &o) noexcept
    {
        return static_cast<CurviCoord &>(Scalar::operator-=(o));
    }
    constexpr decltype(auto) operator*=(CurviCoord const &o) noexcept
    {
        return static_cast<CurviCoord &>(Scalar::operator*=(o));
    }
    constexpr decltype(auto) operator/=(CurviCoord const &o) noexcept
    {
        return static_cast<CurviCoord &>(Scalar::operator/=(o));
    }

    // unary operations
    //
    [[nodiscard]] friend constexpr CurviCoord const &operator+(CurviCoord const &s) noexcept
    {
        return s;
    }
    [[nodiscard]] friend constexpr CurviCoord operator-(CurviCoord const &s) noexcept
    {
        return -Real{ s };
    }

    // binary operations
    //
    [[nodiscard]] friend constexpr CurviCoord operator+(CurviCoord a, CurviCoord const &b) noexcept
    {
        return a += b;
    }
    [[nodiscard]] friend constexpr CurviCoord operator-(CurviCoord a, CurviCoord const &b) noexcept
    {
        return a -= b;
    }
    [[nodiscard]] friend constexpr CurviCoord operator*(CurviCoord a, CurviCoord const &b) noexcept
    {
        return a *= b;
    }
    [[nodiscard]] friend constexpr CurviCoord operator/(CurviCoord a, CurviCoord const &b) noexcept
    {
        return a /= b;
    }

    // pretty print
    //
    template <class CharT, class Traits>
    friend decltype(auto) operator<<(std::basic_ostream<CharT, Traits> &os, CurviCoord const &s)
    {
        return os << static_cast<Scalar const &>(s);
    }
};

static_assert(std::is_standard_layout_v<CurviCoord>);
LIBPIC_END_NAMESPACE
