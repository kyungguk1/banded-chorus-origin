/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#ifndef COMMON_PARTICLE_h
#define COMMON_PARTICLE_h

#include <Utility/Vector.h>
#include <common-config.h>

#include <limits>
#include <ostream>
#include <type_traits>

COMMON_BEGIN_NAMESPACE
/// single particle description
///
struct Particle {
    using Real                      = double;
    static constexpr Real quiet_nan = std::numeric_limits<Real>::quiet_NaN();

    // for delta-f
    struct Delta {
        Real f{ quiet_nan }; // f(0, x(0), v(0))
        Real w{ quiet_nan }; // f(0, x(0), v(0))/g(0, x(0), v(0)) - f_0(x(t), v(t))/g(0, x(0), v(0))
        static constexpr Real fOg{ 1 }; // f(0, x(0), v(0))/g(0, x(0), v(0)),
                                        // where g is the marker particle distribution
    };

    Vector vel{ quiet_nan };   //!< 3-component velocity vector
    Real   pos_x{ quiet_nan }; //!< x-component of position
    Delta  delta{};

    Particle() noexcept = default;
    Particle(Vector const &vel, Real pos_x) noexcept : vel{ vel }, pos_x{ pos_x } {}

    // pretty print
    //
    template <class CharT, class Traits>
    friend decltype(auto) operator<<(std::basic_ostream<CharT, Traits> &os, Particle const &ptl)
    {
        return os << '{' << ptl.vel << ", " << '{' << ptl.pos_x << '}' << ", " << '{' << ptl.delta.f
                  << ", " << ptl.delta.w << '}' << '}';
    }
};

static_assert(sizeof(Particle) == 8 * sizeof(Particle::Real));
static_assert(alignof(Particle) == alignof(Vector));
COMMON_END_NAMESPACE

#endif /* COMMON_PARTICLE_h */
