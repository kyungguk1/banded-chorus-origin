/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#ifndef COMMON_TENSOR_h
#define COMMON_TENSOR_h

#include <Utility/Vector.h>
#include <common-config.h>

#include <ostream>
#include <type_traits>

COMMON_BEGIN_NAMESPACE
/// compact symmetric rank-2 tensor
///
struct alignas(Vector) Tensor {
    using Real     = double;
    using _dummy_t = std::aligned_storage_t<sizeof(Real), alignof(Real)>;

    // tensor elements
    //
    Real     xx{}, yy{}, zz{}; // diagonal
    _dummy_t _dummy1{};
    Real     xy{}, yz{}, zx{}; // off-diag
    _dummy_t _dummy2{};

    // constructors
    //
    constexpr Tensor() noexcept = default;
    constexpr explicit Tensor(Real const v) noexcept : Tensor{ v, v, v, v, v, v } {}
    constexpr Tensor(Real xx, Real yy, Real zz, Real xy, Real yz, Real zx) noexcept
    : xx{ xx }, yy{ yy }, zz{ zz }, xy{ xy }, yz{ yz }, zx{ zx }
    {
    }

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
        return { A.xx * b.x + A.xy * b.y + A.zx * b.z, A.xy * b.x + A.yy * b.y + A.yz * b.z,
                 A.zx * b.x + A.yz * b.y + A.zz * b.z };
    }
    [[nodiscard]] friend constexpr Vector dot(Vector const &b, Tensor const &A) noexcept
    {
        return dot(A, b); // because A^T == A
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
        return Tensor{ a } - b;
    }
    [[nodiscard]] friend constexpr Tensor operator*(Real const &b, Tensor const &a) noexcept
    {
        return a * b;
    }
    [[nodiscard]] friend constexpr Tensor operator/(Real const &a, Tensor const &b) noexcept
    {
        return Tensor{ a } / b;
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
static_assert(alignof(Tensor) == alignof(Vector) && sizeof(Tensor) == 2 * sizeof(Vector),
              "incompatible memory layout");
COMMON_END_NAMESPACE

#endif /* COMMON_TENSOR_h */
