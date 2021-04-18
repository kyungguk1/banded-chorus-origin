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

#ifndef Scalar_h
#define Scalar_h

#include "../Macros.h"
#include "../Predefined.h"

#include <ostream>

HYBRID1D_BEGIN_NAMESPACE
class Scalar {
    Real v{};

public:
    // value access
    //
    constexpr explicit operator Real() const noexcept { return v; }

    // constructors
    //
    constexpr explicit Scalar() noexcept = default;
    constexpr Scalar(Real const v) noexcept : v{v} {}

    // compound operations
    //
    constexpr Scalar &operator+=(Scalar const &o) noexcept
    {
        v += Real{o};
        return *this;
    }
    constexpr Scalar &operator-=(Scalar const &o) noexcept
    {
        v -= Real{o};
        return *this;
    }
    constexpr Scalar &operator*=(Scalar const &o) noexcept
    {
        v *= Real{o};
        return *this;
    }
    constexpr Scalar &operator/=(Scalar const &o) noexcept
    {
        v /= Real{o};
        return *this;
    }

    // unary operations
    //
    [[nodiscard]] friend constexpr Scalar const &operator+(Scalar const &s) noexcept { return s; }
    [[nodiscard]] friend constexpr Scalar operator-(Scalar const &s) noexcept { return -Real{s}; }

    // binary operations
    //
    [[nodiscard]] friend constexpr Scalar operator+(Scalar a, Scalar const &b) noexcept
    {
        return a += b;
    }
    [[nodiscard]] friend constexpr Scalar operator-(Scalar a, Scalar const &b) noexcept
    {
        return a -= b;
    }
    [[nodiscard]] friend constexpr Scalar operator*(Scalar a, Scalar const &b) noexcept
    {
        return a *= b;
    }
    [[nodiscard]] friend constexpr Scalar operator/(Scalar a, Scalar const &b) noexcept
    {
        return a /= b;
    }

    // pretty print
    //
    template <class CharT, class Traits>
    friend decltype(auto) operator<<(std::basic_ostream<CharT, Traits> &os, Scalar const &s)
    {
        return os << Real{s};
    }
};
HYBRID1D_END_NAMESPACE

#endif /* Scalar_h */
