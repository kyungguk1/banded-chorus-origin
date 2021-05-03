/*
 * Copyright (c) 2019, Kyungguk Min
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * - Redistributions of source code must retain the above copyright notice, this
 *  list of conditions and the following disclaimer.
 *
 * - Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef Tensor_h
#define Tensor_h

#include "../Macros.h"
#include "../Predefined.h"
#include "./Vector.h"

#include <ostream>

HYBRID1D_BEGIN_NAMESPACE
/// compact symmetric rank-2 tensor
///
struct Tensor {
    // tensor elements
    //
    Real xx{}, yy{}, zz{}; // diagonal
    Real xy{}, yz{}, zx{}; // off-diag

    // constructors
    //
    constexpr explicit Tensor() noexcept = default;
    constexpr explicit Tensor(Real const v) noexcept
    : xx{ v }, yy{ v }, zz{ v }, xy{ v }, yz{ v }, zx{ v }
    {
    }
    constexpr Tensor(Real xx, Real yy, Real zz, Real xy, Real yz, Real zx) noexcept
    : xx{ xx }, yy{ yy }, zz{ zz }, xy{ xy }, yz{ yz }, zx{ zx }
    {
    }

    // access to lower and upper halves
    //
    [[nodiscard]] Vector &      lo() noexcept { return *reinterpret_cast<Vector *>(&xx); }
    [[nodiscard]] Vector const &lo() const noexcept
    {
        return *reinterpret_cast<Vector const *>(&xx);
    }

    [[nodiscard]] Vector &      hi() noexcept { return *reinterpret_cast<Vector *>(&xy); }
    [[nodiscard]] Vector const &hi() const noexcept
    {
        return *reinterpret_cast<Vector const *>(&xy);
    }

    // compound operations: tensor @= tensor, where @ is one of +, -, *, and / (element-wise)
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

    // compound operations: tensor @= real, where @ is one of +, -, *, and / (applied to all
    // elements)
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
    [[nodiscard]] friend constexpr Tensor operator-(Tensor v) noexcept { return v *= Real{ -1 }; }

    // binary operations: tensor @ {vector|real}, where @ is one of +, -, *, and /
    //
    template <class B>
    [[nodiscard]] friend constexpr Tensor operator+(Tensor a, B const &b) noexcept
    {
        return a += b;
    }
    template <class B>
    [[nodiscard]] friend constexpr Tensor operator-(Tensor a, B const &b) noexcept
    {
        return a -= b;
    }
    template <class B>
    [[nodiscard]] friend constexpr Tensor operator*(Tensor a, B const &b) noexcept
    {
        return a *= b;
    }
    template <class B>
    [[nodiscard]] friend constexpr Tensor operator/(Tensor a, B const &b) noexcept
    {
        return a /= b;
    }

    // binary operations: real @ tensor, where @ is one of +, -, *, and /
    //
    [[nodiscard]] friend constexpr Tensor operator+(Real const &b, Tensor a) noexcept
    {
        return a += b;
    }
    [[nodiscard]] friend constexpr Tensor operator-(Real a, Tensor const &b) noexcept
    {
        return Tensor{ a } -= b;
    }
    [[nodiscard]] friend constexpr Tensor operator*(Real const &b, Tensor a) noexcept
    {
        return a *= b;
    }
    [[nodiscard]] friend constexpr Tensor operator/(Real a, Tensor const &b) noexcept
    {
        return Tensor{ a } /= b;
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
HYBRID1D_END_NAMESPACE

#endif /* Tensor_h */
