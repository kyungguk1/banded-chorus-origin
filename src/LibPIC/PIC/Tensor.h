/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC/Config.h>
#include <PIC/Vector.h>

#include <functional>
#include <ostream>
#include <type_traits>

LIBPIC_BEGIN_NAMESPACE
/// compact symmetric rank-2 tensor
///
struct alignas(Vector) Tensor {
    using Real = double;

    // tensor elements
    //
    Real xx{}, yy{}, zz{}; // diagonal
    Real xy{}, yz{}, zx{}; // off-diag

    [[nodiscard]] constexpr Real &m11() noexcept { return xx; }
    [[nodiscard]] constexpr Real &m12() noexcept { return xy; }
    [[nodiscard]] constexpr Real &m13() noexcept { return zx; }
    [[nodiscard]] constexpr Real &m21() noexcept { return xy; }
    [[nodiscard]] constexpr Real &m22() noexcept { return yy; }
    [[nodiscard]] constexpr Real &m23() noexcept { return yz; }
    [[nodiscard]] constexpr Real &m31() noexcept { return zx; }
    [[nodiscard]] constexpr Real &m32() noexcept { return yz; }
    [[nodiscard]] constexpr Real &m33() noexcept { return zz; }

    [[nodiscard]] constexpr Real const &m11() const noexcept { return xx; }
    [[nodiscard]] constexpr Real const &m12() const noexcept { return xy; }
    [[nodiscard]] constexpr Real const &m13() const noexcept { return zx; }
    [[nodiscard]] constexpr Real const &m21() const noexcept { return xy; }
    [[nodiscard]] constexpr Real const &m22() const noexcept { return yy; }
    [[nodiscard]] constexpr Real const &m23() const noexcept { return yz; }
    [[nodiscard]] constexpr Real const &m31() const noexcept { return zx; }
    [[nodiscard]] constexpr Real const &m32() const noexcept { return yz; }
    [[nodiscard]] constexpr Real const &m33() const noexcept { return zz; }

    // constructors
    //
    constexpr Tensor() noexcept = default;
    constexpr explicit Tensor(Real const v) noexcept
    : Tensor{ v, v, v, v, v, v } {}
    constexpr Tensor(Real xx, Real yy, Real zz, Real xy, Real yz, Real zx) noexcept
    : xx{ xx }, yy{ yy }, zz{ zz }, xy{ xy }, yz{ yz }, zx{ zx } {}

    [[nodiscard]] static constexpr auto identity() noexcept { return Tensor{ 1, 1, 1, 0, 0, 0 }; }

    // access to lower and upper halves as a vector
    //
    [[nodiscard]] Vector       &lo() noexcept { return *reinterpret_cast<Vector *>(&xx); }
    [[nodiscard]] Vector const &lo() const noexcept
    {
        return *reinterpret_cast<Vector const *>(&xx);
    }

    [[nodiscard]] Vector       &hi() noexcept { return *reinterpret_cast<Vector *>(&xy); }
    [[nodiscard]] Vector const &hi() const noexcept
    {
        return *reinterpret_cast<Vector const *>(&xy);
    }

    // tensor calculus
    //
    [[nodiscard]] friend constexpr Vector dot(Tensor const &A, Vector const &b) noexcept
    {
        return {
            A.xx * b.x + A.xy * b.y + A.zx * b.z,
            A.xy * b.x + A.yy * b.y + A.yz * b.z,
            A.zx * b.x + A.yz * b.y + A.zz * b.z,
        };
    }
    [[nodiscard]] friend constexpr Vector dot(Vector const &b, Tensor const &A) noexcept
    {
        return dot(A, b); // because A^T == A
    }
    [[nodiscard]] friend constexpr Real trace(Tensor const &A) noexcept
    {
        return A.xx + A.yy + A.zz;
    }
    [[nodiscard]] friend constexpr Tensor const &transpose(Tensor const &A) noexcept
    {
        return A;
    }
    [[nodiscard]] friend constexpr Real det(Tensor const &A) noexcept
    {
        return (A.xx * A.yy * A.zz + 2 * A.xy * A.yz * A.zx)
             - (A.xx * A.yz * A.yz + A.yy * A.zx * A.zx + A.xy * A.xy * A.zz);
    }
    [[nodiscard]] friend constexpr Tensor inv(Tensor const &A) noexcept
    {
        Tensor inv{
            A.yy * A.zz - A.yz * A.yz, A.xx * A.zz - A.zx * A.zx, A.xx * A.yy - A.xy * A.xy,
            A.yz * A.zx - A.xy * A.zz, A.xy * A.zx - A.xx * A.yz, A.xy * A.yz - A.yy * A.zx
        };
        inv /= det(A);
        return inv;
    }

    // left-fold: applies to all elements
    // the signature of BinaryOp is "Init(Init, Real)"
    //
    template <class Init, class BinaryOp,
              std::enable_if_t<std::is_invocable_r_v<Init, BinaryOp, Init, Real>, int> = 0>
    [[nodiscard]] constexpr auto fold(Init init, BinaryOp &&f) const
        noexcept(noexcept(std::invoke(f, init, Real{})))
    {
        return f(f(f(f(f(f(init, xx), yy), zz), xy), yz), zx);
    }

    // compound operations: tensor @= tensor, where @ is one of +, -, *, and /
    // operation is element-wise
    //
    constexpr Tensor &operator+=(Tensor const &v) noexcept
    {
        xx += v.xx;
        yy += v.yy;
        zz += v.zz;
        xy += v.xy;
        yz += v.yz;
        zx += v.zx;
        return *this;
    }
    constexpr Tensor &operator-=(Tensor const &v) noexcept
    {
        xx -= v.xx;
        yy -= v.yy;
        zz -= v.zz;
        xy -= v.xy;
        yz -= v.yz;
        zx -= v.zx;
        return *this;
    }
    constexpr Tensor &operator*=(Tensor const &v) noexcept
    {
        xx *= v.xx;
        yy *= v.yy;
        zz *= v.zz;
        xy *= v.xy;
        yz *= v.yz;
        zx *= v.zx;
        return *this;
    }
    constexpr Tensor &operator/=(Tensor const &v) noexcept
    {
        xx /= v.xx;
        yy /= v.yy;
        zz /= v.zz;
        xy /= v.xy;
        yz /= v.yz;
        zx /= v.zx;
        return *this;
    }

    // scalar-tensor compound operations: tensor @= real, where @ is one of +, -, *, and /
    // operation with scalar is distributed to all elements
    //
    constexpr Tensor &operator+=(Real const &s) noexcept
    {
        xx += s;
        yy += s;
        zz += s;
        xy += s;
        yz += s;
        zx += s;
        return *this;
    }
    constexpr Tensor &operator-=(Real const &s) noexcept
    {
        xx -= s;
        yy -= s;
        zz -= s;
        xy -= s;
        yz -= s;
        zx -= s;
        return *this;
    }
    constexpr Tensor &operator*=(Real const &s) noexcept
    {
        xx *= s;
        yy *= s;
        zz *= s;
        xy *= s;
        yz *= s;
        zx *= s;
        return *this;
    }
    constexpr Tensor &operator/=(Real const &s) noexcept
    {
        xx /= s;
        yy /= s;
        zz /= s;
        xy /= s;
        yz /= s;
        zx /= s;
        return *this;
    }

    // unary operations
    //
    [[nodiscard]] friend constexpr Tensor const &operator+(Tensor const &v) noexcept { return v; }
    [[nodiscard]] friend constexpr Tensor        operator-(Tensor v) noexcept
    {
        v *= Real{ -1 };
        return v;
    }

    // binary operations: tensor @ {tensor|real}, where @ is one of +, -, *, and /
    //
    template <class B>
    [[nodiscard]] friend constexpr Tensor operator+(Tensor a, B const &b) noexcept
    {
        a += b;
        return a;
    }
    template <class B>
    [[nodiscard]] friend constexpr Tensor operator-(Tensor a, B const &b) noexcept
    {
        a -= b;
        return a;
    }
    template <class B>
    [[nodiscard]] friend constexpr Tensor operator*(Tensor a, B const &b) noexcept
    {
        a *= b;
        return a;
    }
    template <class B>
    [[nodiscard]] friend constexpr Tensor operator/(Tensor a, B const &b) noexcept
    {
        a /= b;
        return a;
    }

    // binary operations: real @ tensor, where @ is one of +, -, *, and /
    //
    [[nodiscard]] friend constexpr Tensor operator+(Real const &b, Tensor const &a) noexcept
    {
        return a + b;
    }
    [[nodiscard]] friend constexpr Tensor operator-(Real const &a, Tensor const &b) noexcept
    {
        Tensor A{ a };
        A -= b;
        return A;
    }
    [[nodiscard]] friend constexpr Tensor operator*(Real const &b, Tensor const &a) noexcept
    {
        return a * b;
    }
    [[nodiscard]] friend constexpr Tensor operator/(Real const &a, Tensor const &b) noexcept
    {
        Tensor A{ a };
        A /= b;
        return A;
    }

    // pretty print
    //
    template <class CharT, class Traits>
    friend decltype(auto) operator<<(std::basic_ostream<CharT, Traits> &os, Tensor const &v)
    {
        return os << '{' << v.xx << ", " << v.yy << ", " << v.zz << ", " << v.xy << ", " << v.yz
                  << ", " << v.zx << '}';
    }
};

// make sure that memory layout of Tensor and Vector are compatible
//
static_assert(alignof(Tensor) == alignof(Vector));
static_assert(sizeof(Tensor) == 2 * sizeof(Vector));
static_assert(std::is_standard_layout_v<Tensor>);
LIBPIC_END_NAMESPACE
