/*
 * Copyright (c) 2021-2022, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC/Config.h>
#include <PIC/Scalar.h>
#include <PIC/Tensor.h>
#include <PIC/Vector.h>

#include <ostream>

LIBPIC_NAMESPACE_BEGIN(1)
/// Symmetric rank-2 four-tensor
///
struct FourTensor {
    using Real = double;

    // tensor elements
    //
    Scalar tt{}; // time component
    Vector ts{}; // mixed components
    Tensor ss{}; // space components

    // constructors
    //
    constexpr FourTensor() noexcept = default;
    constexpr explicit FourTensor(Real const v) noexcept
    : tt{ v }, ts{ v }, ss{ v } {}
    constexpr FourTensor(Scalar const &tt, Vector const &ts, Tensor const &ss) noexcept
    : tt{ tt }, ts{ ts }, ss{ ss }
    {
    }

    // compound operations: tensor @= tensor, where @ is one of +, -, *, and /
    // operation is element-wise
    //
    constexpr FourTensor &operator+=(FourTensor const &v) noexcept
    {
        tt += v.tt;
        ts += v.ts;
        ss += v.ss;
        return *this;
    }
    constexpr FourTensor &operator-=(FourTensor const &v) noexcept
    {
        tt -= v.tt;
        ts -= v.ts;
        ss -= v.ss;
        return *this;
    }
    constexpr FourTensor &operator*=(FourTensor const &v) noexcept
    {
        tt *= v.tt;
        ts *= v.ts;
        ss *= v.ss;
        return *this;
    }
    constexpr FourTensor &operator/=(FourTensor const &v) noexcept
    {
        tt /= v.tt;
        ts /= v.ts;
        ss /= v.ss;
        return *this;
    }

    // compound operations: tensor @= real, where @ is one of +, -, *, and /
    // operation with scalar is distributed to all elements
    //
    constexpr FourTensor &operator+=(Real const &v) noexcept
    {
        tt += v;
        ts += v;
        ss += v;
        return *this;
    }
    constexpr FourTensor &operator-=(Real const &v) noexcept
    {
        tt -= v;
        ts -= v;
        ss -= v;
        return *this;
    }
    constexpr FourTensor &operator*=(Real const &v) noexcept
    {
        tt *= v;
        ts *= v;
        ss *= v;
        return *this;
    }
    constexpr FourTensor &operator/=(Real const &v) noexcept
    {
        tt /= v;
        ts /= v;
        ss /= v;
        return *this;
    }

    // unary operations
    //
    [[nodiscard]] friend constexpr FourTensor const &operator+(FourTensor const &v) noexcept
    {
        return v;
    }
    [[nodiscard]] friend constexpr FourTensor operator-(FourTensor v) noexcept
    {
        v *= Real{ -1 };
        return v;
    }

    // binary operations: tensor @ {tensor|real}, where @ is one of +, -, *, and /
    //
    template <class B>
    [[nodiscard]] friend constexpr FourTensor operator+(FourTensor a, B const &b) noexcept
    {
        a += b;
        return a;
    }
    template <class B>
    [[nodiscard]] friend constexpr FourTensor operator-(FourTensor a, B const &b) noexcept
    {
        a -= b;
        return a;
    }
    template <class B>
    [[nodiscard]] friend constexpr FourTensor operator*(FourTensor a, B const &b) noexcept
    {
        a *= b;
        return a;
    }
    template <class B>
    [[nodiscard]] friend constexpr FourTensor operator/(FourTensor a, B const &b) noexcept
    {
        a /= b;
        return a;
    }

    // binary operations: real @ tensor, where @ is one of +, -, *, and /
    //
    [[nodiscard]] friend constexpr FourTensor operator+(Real const &b, FourTensor const &a) noexcept
    {
        return a + b;
    }
    [[nodiscard]] friend constexpr FourTensor operator-(Real const &a, FourTensor const &b) noexcept
    {
        FourTensor A{ a };
        A -= b;
        return A;
    }
    [[nodiscard]] friend constexpr FourTensor operator*(Real const &b, FourTensor const &a) noexcept
    {
        return a * b;
    }
    [[nodiscard]] friend constexpr FourTensor operator/(Real const &a, FourTensor const &b) noexcept
    {
        FourTensor A{ a };
        A /= b;
        return A;
    }

    // pretty print
    //
    template <class CharT, class Traits>
    friend decltype(auto) operator<<(std::basic_ostream<CharT, Traits> &os, FourTensor const &v)
    {
        return os << '{' << v.tt << ", " << v.ts << ", " << v.ss << '}';
    }
};
LIBPIC_NAMESPACE_END(1)
