/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC/Config.h>

#include <ostream>
#include <type_traits>

LIBPIC_BEGIN_NAMESPACE
struct alignas(32) Vector {
    using Real    = double;
    using Padding = std::aligned_storage_t<sizeof(Real), alignof(Real)>;

    // vector elements
    //
    Real    x{};
    Real    y{};
    Real    z{};
    Padding _{};

    // constructors
    //
    constexpr Vector() noexcept = default;
    constexpr explicit Vector(Real const v) noexcept : Vector{ v, v, v } {}
    constexpr Vector(Real const x, Real const y, Real const z) noexcept : x{ x }, y{ y }, z{ z } {}

    // vector calculus
    //
    [[nodiscard]] friend constexpr Real dot(Vector const &A, Vector const &B) noexcept
    {
        return A.x * B.x + A.y * B.y + A.z * B.z;
    }
    [[nodiscard]] friend constexpr Vector cross(Vector const &A, Vector const &B) noexcept
    {
        return { A.y * B.z - A.z * B.y, A.z * B.x - A.x * B.z, A.x * B.y - A.y * B.x };
    }

    // left-fold: applies to all elements
    // the signature of BinaryOp is "Init(Init, Real)"
    //
    template <class Init, class BinaryOp,
              std::enable_if_t<std::is_invocable_r_v<Init, BinaryOp, Init, Real>, int> = 0>
    [[nodiscard]] constexpr auto fold(Init init, BinaryOp &&f) const
        noexcept(noexcept(std::invoke(f, init, Real{})))
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

    // scalar-vector compound operations: vector @= real, where @ is one of +, -, *, and /
    // operation with scalar is distributed to all elements
    //
    constexpr Vector &operator+=(Real const &s) noexcept
    {
        x += s;
        y += s;
        z += s;
        return *this;
    }
    constexpr Vector &operator-=(Real const &s) noexcept
    {
        x -= s;
        y -= s;
        z -= s;
        return *this;
    }
    constexpr Vector &operator*=(Real const &s) noexcept
    {
        x *= s;
        y *= s;
        z *= s;
        return *this;
    }
    constexpr Vector &operator/=(Real const &s) noexcept
    {
        x /= s;
        y /= s;
        z /= s;
        return *this;
    }

    // unary operations
    //
    [[nodiscard]] friend constexpr Vector const &operator+(Vector const &v) noexcept { return v; }
    [[nodiscard]] friend constexpr Vector        operator-(Vector v) noexcept
    {
        v *= Real{ -1 };
        return v;
    }

    // binary operations: vector @ {vector|real}, where @ is one of +, -, *, and /
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

    // binary operations: real @ vector, where @ is one of +, -, *, and /
    //
    [[nodiscard]] friend constexpr Vector operator+(Real const &b, Vector const &a) noexcept
    {
        return a + b;
    }
    [[nodiscard]] friend constexpr Vector operator-(Real const &a, Vector const &b) noexcept
    {
        return Vector{ a } - b;
    }
    [[nodiscard]] friend constexpr Vector operator*(Real const &b, Vector const &a) noexcept
    {
        return a * b;
    }
    [[nodiscard]] friend constexpr Vector operator/(Real const &a, Vector const &b) noexcept
    {
        return Vector{ a } / b;
    }

    // pretty print
    //
    template <class CharT, class Traits>
    friend decltype(auto) operator<<(std::basic_ostream<CharT, Traits> &os, Vector const &v)
    {
        return os << '{' << v.x << ", " << v.y << ", " << v.z << '}';
    }
};

static_assert(32 == sizeof(Vector));
static_assert(32 == alignof(Vector));
static_assert(std::is_standard_layout_v<Vector>);
LIBPIC_END_NAMESPACE
