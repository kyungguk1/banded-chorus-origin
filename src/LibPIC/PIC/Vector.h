/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC/Config.h>

#include <complex>
#include <functional>
#include <ostream>
#include <type_traits>

LIBPIC_BEGIN_NAMESPACE
/// Generic Vector
template <class Type = double>
struct GenericVector {
    using Vector     = GenericVector<Type>;
    using value_type = Type;

    // vector elements
    //
    Type x{};
    Type y{};
    Type z{};

    // constructors
    //
    constexpr GenericVector() noexcept = default;
    constexpr explicit GenericVector(Type const &v) noexcept
    : GenericVector{ v, v, v } {}
    constexpr GenericVector(Type const &x, Type const &y, Type const &z) noexcept
    : x{ x }, y{ y }, z{ z } {}

    // vector calculus
    //
    [[nodiscard]] friend constexpr Type dot(Vector const &A, Vector const &B) noexcept
    {
        return A.x * B.x + A.y * B.y + A.z * B.z;
    }
    [[nodiscard]] friend constexpr Vector cross(Vector const &A, Vector const &B) noexcept
    {
        return { A.y * B.z - A.z * B.y, A.z * B.x - A.x * B.z, A.x * B.y - A.y * B.x };
    }

    // left-fold: applies to all elements
    // the signature of BinaryOp is "Init(Init, Type)"
    //
    template <class Init, class BinaryOp,
              std::enable_if_t<std::is_invocable_r_v<Init, BinaryOp, Init, Type>, int> = 0>
    [[nodiscard]] constexpr auto fold(Init init, BinaryOp &&f) const
        noexcept(std::is_nothrow_invocable_r_v<Init, BinaryOp, Init, Type>)
    {
        return f(f(f(init, x), y), z);
    }

    // compound operations: vector @= vector, where @ is one of +, -, *, and /
    // operation is element-wise
    //
    constexpr Vector &operator+=(Vector const &v) noexcept
    {
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }
    constexpr Vector &operator-=(Vector const &v) noexcept
    {
        x -= v.x;
        y -= v.y;
        z -= v.z;
        return *this;
    }
    constexpr Vector &operator*=(Vector const &v) noexcept
    {
        x *= v.x;
        y *= v.y;
        z *= v.z;
        return *this;
    }
    constexpr Vector &operator/=(Vector const &v) noexcept
    {
        x /= v.x;
        y /= v.y;
        z /= v.z;
        return *this;
    }

    // scalar-vector compound operations: vector @= type, where @ is one of +, -, *, and /
    // operation with scalar is distributed to all elements
    //
    constexpr Vector &operator+=(Type const &s) noexcept
    {
        x += s;
        y += s;
        z += s;
        return *this;
    }
    constexpr Vector &operator-=(Type const &s) noexcept
    {
        x -= s;
        y -= s;
        z -= s;
        return *this;
    }
    constexpr Vector &operator*=(Type const &s) noexcept
    {
        x *= s;
        y *= s;
        z *= s;
        return *this;
    }
    constexpr Vector &operator/=(Type const &s) noexcept
    {
        x /= s;
        y /= s;
        z /= s;
        return *this;
    }

    // unary operations
    //
    [[nodiscard]] friend constexpr Vector const &operator+(Vector const &v) noexcept { return v; }
    [[nodiscard]] friend constexpr Vector        operator-(Vector const &v) noexcept { return Vector{} - v; }

    // binary operations: vector @ {vector|type}, where @ is one of +, -, *, and /
    //
    template <class B>
    [[nodiscard]] friend constexpr Vector operator+(Vector a, B const &b) noexcept
    {
        a += b;
        return a;
    }
    template <class B>
    [[nodiscard]] friend constexpr Vector operator-(Vector a, B const &b) noexcept
    {
        a -= b;
        return a;
    }
    template <class B>
    [[nodiscard]] friend constexpr Vector operator*(Vector a, B const &b) noexcept
    {
        a *= b;
        return a;
    }
    template <class B>
    [[nodiscard]] friend constexpr Vector operator/(Vector a, B const &b) noexcept
    {
        a /= b;
        return a;
    }

    // binary operations: type @ vector, where @ is one of +, -, *, and /
    //
    [[nodiscard]] friend constexpr Vector operator+(Type const &b, Vector const &a) noexcept
    {
        return a + b;
    }
    [[nodiscard]] friend constexpr Vector operator-(Type const &a, Vector const &b) noexcept
    {
        Vector A{ a };
        A -= b;
        return A;
    }
    [[nodiscard]] friend constexpr Vector operator*(Type const &b, Vector const &a) noexcept
    {
        return a * b;
    }
    [[nodiscard]] friend constexpr Vector operator/(Type const &a, Vector const &b) noexcept
    {
        Vector A{ a };
        A /= b;
        return A;
    }

    // pretty print
    //
    template <class CharT, class Traits>
    friend decltype(auto) operator<<(std::basic_ostream<CharT, Traits> &os, Vector const &v)
    {
        return os << '{' << v.x << ", " << v.y << ", " << v.z << '}';
    }
};

// real vector
using Vector = GenericVector<double>;
static_assert(24 == sizeof(Vector));
static_assert(8 == alignof(Vector));
static_assert(std::is_standard_layout_v<Vector>);

// complex vector
using ComplexVector = GenericVector<std::complex<double>>;
using namespace std::literals::complex_literals;
LIBPIC_END_NAMESPACE
