/*
 * Copyright (c) 2020, Kyungguk Min
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

#ifndef Range_h
#define Range_h

#include "../Macros.h"
#include "../Predefined.h"

#include <ostream>

HYBRID1D_BEGIN_NAMESPACE
/// represents a range between two points, a and b.
///
struct [[nodiscard]] Range {
    Real loc; //!< beginning of the range.
    Real len; //!< length of the interval; must be non-negative.
    //
    [[nodiscard]] constexpr Real min() const noexcept { return loc; }
    [[nodiscard]] constexpr Real max() const noexcept { return loc + len; }

    /// return true if a point, x, is contained in [a, b)
    ///
    [[nodiscard]] constexpr bool is_member(Real const x) const noexcept
    {
        return x >= min() && x < max();
    }

    // compound operations
    //
    constexpr Range &operator+=(Real const &o) noexcept
    {
        loc += o;
        return *this;
    }
    constexpr Range &operator-=(Real const &o) noexcept { return *this += -o; }
    constexpr Range &operator*=(Real const &o) noexcept
    {
        if (o >= 0) {
            loc *= o;
            len *= o;
        } else {
            loc = o * max();
            len *= -o;
        }
        return *this;
    }
    constexpr Range &operator/=(Real const &o) noexcept { return *this *= 1 / o; }

    // unary operations
    //
    [[nodiscard]] friend constexpr Range const &operator+(Range const &s) noexcept { return s; }
    [[nodiscard]] friend constexpr Range        operator-(Range s) noexcept
    {
        s *= -1;
        return s;
    }

    // binary operations
    //
    [[nodiscard]] friend constexpr Range operator+(Range a, Real const &b) noexcept
    {
        a += b;
        return a;
    }
    [[nodiscard]] friend constexpr Range operator-(Range a, Real const &b) noexcept
    {
        a -= b;
        return a;
    }
    [[nodiscard]] friend constexpr Range operator*(Range a, Real const &b) noexcept
    {
        a *= b;
        return a;
    }
    [[nodiscard]] friend constexpr Range operator/(Range a, Real const &b) noexcept
    {
        a /= b;
        return a;
    }

    [[nodiscard]] friend constexpr Range operator+(Real const &b, Range a) noexcept
    {
        a += b;
        return a;
    }
    [[nodiscard]] friend constexpr Range operator-(Real const &b, Range const &a) noexcept
    {
        return (-a) + b;
    }
    [[nodiscard]] friend constexpr Range operator*(Real const &b, Range a) noexcept
    {
        a *= b;
        return a;
    }

private:
    // pretty print
    //
    template <class CharT, class Traits>
    friend decltype(auto) operator<<(std::basic_ostream<CharT, Traits> &os, Range const &r)
    {
        return os << '{' << r.min() << ", " << r.max() << '}';
    }
};
HYBRID1D_END_NAMESPACE

#endif /* Range_h */
