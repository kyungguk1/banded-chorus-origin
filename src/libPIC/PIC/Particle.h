/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC/Config.h>
#include <PIC/Vector.h>

#include <array>
#include <limits>
#include <ostream>
#include <type_traits>

LIBPIC_BEGIN_NAMESPACE
/// single particle description
///
struct Particle {
    using Real                      = double;
    using pad_t                     = std::array<unsigned char, sizeof(Real)>;
    static constexpr Real quiet_nan = std::numeric_limits<Real>::quiet_NaN();

    // for delta-f
    struct PSD {
        Real full_f{ quiet_nan }; // f(0, x(0), v(0))
        Real weight{ quiet_nan }; // f(0, x(0), v(0))/g(0, x(0), v(0)) - f_0(x(t), v(t))/g(0, x(0), v(0))

        // f(0, x(0), v(0))/g(0, x(0), v(0)), where g is the marker particle distribution
        static constexpr Real fOg{ 1 };
    };

    Vector vel{ quiet_nan };   //!< 3-component velocity vector
    Real   pos_x{ quiet_nan }; //!< x-component of position
    PSD    psd{};
    long   id{ -1 }; //!< particle identifier
    pad_t  padding{};

    Particle() noexcept = default;
    Particle(Vector const &vel, Real pos_x) noexcept
    : vel{ vel }, pos_x{ pos_x }, id{ next_id() }
    {
    }

private:
    [[nodiscard]] static long next_id() noexcept
    {
        thread_local static long next_id{ 0 };
        return next_id++;
    }

    // pretty print
    //
    template <class CharT, class Traits>
    friend decltype(auto) operator<<(std::basic_ostream<CharT, Traits> &os, Particle const &ptl)
    {
        return os << '{' << ptl.vel << ", " << '{' << ptl.pos_x << '}' << ", "
                  << '{' << ptl.psd.full_f << ", " << ptl.psd.weight << '}' << '}';
    }
};
static_assert(sizeof(Particle) == 8 * sizeof(Particle::Real));
static_assert(alignof(Particle) == alignof(Vector));
static_assert(std::is_standard_layout_v<Particle>);
LIBPIC_END_NAMESPACE
