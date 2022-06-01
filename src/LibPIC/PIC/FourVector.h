/*
 * Copyright (c) 2021-2022, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC/Config.h>
#include <PIC/Scalar.h>
#include <PIC/Vector.h>

#include <ostream>
#include <type_traits>

LIBPIC_NAMESPACE_BEGIN(1)
/// Four vector
///
struct FourVector {
    using Real = double;

    // vector elements
    //
    Scalar t{}; //!< time component
    Vector s{}; //!< space components

    // constructors
    //
    constexpr FourVector() noexcept = default;
    constexpr explicit FourVector(Real v) noexcept
    : t{ v }, s{ v } {}
    constexpr FourVector(Scalar const &t, Vector const &s) noexcept
    : t{ t }, s{ s } {}

    // compound operations: vector @= vector, where @ is one of +, -, *, and /
    // operation is element-wise
    //
    constexpr FourVector &operator+=(FourVector const &v) noexcept
    {
        t += v.t;
        s += v.s;
        return *this;
    }
    constexpr FourVector &operator-=(FourVector const &v) noexcept
    {
        t -= v.t;
        s -= v.s;
        return *this;
    }
    constexpr FourVector &operator*=(FourVector const &v) noexcept
    {
        t *= v.t;
        s *= v.s;
        return *this;
    }
    constexpr FourVector &operator/=(FourVector const &v) noexcept
    {
        t /= v.t;
        s /= v.s;
        return *this;
    }

    // compound operations: vector @= real, where @ is one of +, -, *, and /
    // operation with scalar is distributed to all elements
    //
    constexpr FourVector &operator+=(Real const &v) noexcept
    {
        t += v;
        s += v;
        return *this;
    }
    constexpr FourVector &operator-=(Real const &v) noexcept
    {
        t -= v;
        s -= v;
        return *this;
    }
    constexpr FourVector &operator*=(Real const &v) noexcept
    {
        t *= v;
        s *= v;
        return *this;
    }
    constexpr FourVector &operator/=(Real const &v) noexcept
    {
        t /= v;
        s /= v;
        return *this;
    }

    // unary operations
    //
    [[nodiscard]] friend constexpr FourVector const &operator+(FourVector const &v) noexcept
    {
        return v;
    }
    [[nodiscard]] friend constexpr FourVector operator-(FourVector v) noexcept
    {
        v *= Real{ -1 };
        return v;
    }

    // binary operations: vector @ {vector|real}, where @ is one of +, -, *, and /
    //
    template <class B>
    [[nodiscard]] friend constexpr FourVector operator+(FourVector a, B const &b) noexcept
    {
        a += b;
        return a;
    }
    template <class B>
    [[nodiscard]] friend constexpr FourVector operator-(FourVector a, B const &b) noexcept
    {
        a -= b;
        return a;
    }
    template <class B>
    [[nodiscard]] friend constexpr FourVector operator*(FourVector a, B const &b) noexcept
    {
        a *= b;
        return a;
    }
    template <class B>
    [[nodiscard]] friend constexpr FourVector operator/(FourVector a, B const &b) noexcept
    {
        a /= b;
        return a;
    }

    // binary operations: real @ vector, where @ is one of +, -, *, and /
    //
    [[nodiscard]] friend constexpr FourVector operator+(Real const &b, FourVector const &a) noexcept
    {
        return a + b;
    }
    [[nodiscard]] friend constexpr FourVector operator-(Real const &a, FourVector const &b) noexcept
    {
        FourVector A{ a };
        A -= b;
        return A;
    }
    [[nodiscard]] friend constexpr FourVector operator*(Real const &b, FourVector const &a) noexcept
    {
        return a * b;
    }
    [[nodiscard]] friend constexpr FourVector operator/(Real const &a, FourVector const &b) noexcept
    {
        FourVector A{ a };
        A /= b;
        return A;
    }

    // pretty print
    //
    template <class CharT, class Traits>
    friend decltype(auto) operator<<(std::basic_ostream<CharT, Traits> &os, FourVector const &v)
    {
        return os << '{' << v.t << ", " << v.s << '}';
    }
};
LIBPIC_NAMESPACE_END(1)
