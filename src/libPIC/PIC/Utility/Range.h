/*
 * Copyright (c) 2020-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC-config.h>
#include <Predefined.h>

#include <ostream>

LIBPIC_BEGIN_NAMESPACE
/// represents a range between two points, a and b.
///
struct alignas(16) Range {
    Real loc; //!< beginning of the range.
    Real len; //!< length of the interval; must be non-negative.

    /// min and max of the range
    ///
    [[nodiscard]] constexpr Real min() const noexcept { return loc; }
    [[nodiscard]] constexpr Real max() const noexcept { return loc + len; }
    [[nodiscard]] constexpr Real mean() const noexcept { return loc + len / 2; }

    /// return true if a point, x, is contained in [a, b)
    ///
    [[nodiscard]] constexpr bool is_member(Real const x) const noexcept
    {
        return x >= min() && x < max();
    }

    // compound operations
    //
    constexpr Range &operator-=(Real const &o) noexcept { return *this += -o; }
    constexpr Range &operator/=(Real const &o) noexcept { return *this *= 1 / o; }
    constexpr Range &operator+=(Real const &o) noexcept
    {
        loc += o;
        return *this;
    }
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

static_assert(16 == sizeof(Range));
static_assert(16 == alignof(Range));
static_assert(std::is_standard_layout_v<Range>);
LIBPIC_END_NAMESPACE
