//
// Copyright (c) 2019, Kyungguk Min
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// - Redistributions of source code must retain the above copyright notice, this
//  list of conditions and the following disclaimer.
//
// - Redistributions in binary form must reproduce the above copyright notice,
//  this list of conditions and the following disclaimer in the documentation
//  and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#ifndef Vector_h
#define Vector_h

#include "../Macros.h"
#include "../Predefined.h"

#include <ostream>

HYBRID1D_BEGIN_NAMESPACE
struct Vector {
    // vector elements
    //
    Real x{};
    Real y{};
    Real z{};

    // constructors
    //
    constexpr explicit Vector() noexcept = default;
    constexpr explicit Vector(Real const v) noexcept : x{v}, y{v}, z{v} {}
    constexpr Vector(Real const x, Real const y, Real const z) noexcept : x{x}, y{y}, z{z} {}

    // vector calculus
    //
    [[nodiscard]] friend constexpr Real dot(Vector const &A, Vector const &B) noexcept
    {
        return A.x * B.x + A.y * B.y + A.z * B.z;
    }
    [[nodiscard]] friend constexpr Vector cross(Vector const &A, Vector const &B) noexcept
    {
        return {A.y * B.z - A.z * B.y, A.z * B.x - A.x * B.z, A.x * B.y - A.y * B.x};
    }

    // compound operations: vector @= vector, where @ is one of +, -, *, and / (element-wise)
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

    // compound operations: vector @= real, where @ is one of +, -, *, and / (applied to all
    // elements)
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
        v *= Real{-1};
        return v; // NRVO
    }

    // binary operations: vector @ {vector|real}, where @ is one of +, -, *, and /
    //
    template <class B>
    [[nodiscard]] friend constexpr Vector operator+(Vector a, B const &b) noexcept
    {
        return a += b;
    }
    template <class B>
    [[nodiscard]] friend constexpr Vector operator-(Vector a, B const &b) noexcept
    {
        return a -= b;
    }
    template <class B>
    [[nodiscard]] friend constexpr Vector operator*(Vector a, B const &b) noexcept
    {
        return a *= b;
    }
    template <class B>
    [[nodiscard]] friend constexpr Vector operator/(Vector a, B const &b) noexcept
    {
        return a /= b;
    }

    // binary operations: real @ vector, where @ is one of +, -, *, and /
    //
    [[nodiscard]] friend constexpr Vector operator+(Real const &b, Vector a) noexcept
    {
        return a += b;
    }
    [[nodiscard]] friend constexpr Vector operator-(Real a, Vector const &b) noexcept
    {
        return Vector{a} -= b;
    }
    [[nodiscard]] friend constexpr Vector operator*(Real const &b, Vector a) noexcept
    {
        return a *= b;
    }
    [[nodiscard]] friend constexpr Vector operator/(Real a, Vector const &b) noexcept
    {
        return Vector{a} /= b;
    }

    // pretty print
    //
    template <class CharT, class Traits>
    friend decltype(auto) operator<<(std::basic_ostream<CharT, Traits> &os, Vector const &v)
    {
        return os << '{' << v.x << ", " << v.y << ", " << v.z << '}';
    }
};
HYBRID1D_END_NAMESPACE

#endif /* Vector_h */
