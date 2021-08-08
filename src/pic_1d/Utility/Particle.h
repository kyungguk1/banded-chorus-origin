/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#ifndef Particle_h
#define Particle_h

#include "../Macros.h"
#include "../Predefined.h"
#include "./Vector.h"

#include <limits>
#include <ostream>

PIC1D_BEGIN_NAMESPACE
/// single particle description
///
struct Particle {
    static constexpr Real quiet_nan = std::numeric_limits<Real>::quiet_NaN();

    Vector vel{ quiet_nan };   //!< 3-component velocity vector
    Real   pos_x{ quiet_nan }; //!< x-component of position

    // for delta-f
    //
    Real f{ quiet_nan }; // f(0, x(0), v(0))
    Real w{ quiet_nan }; // f(0, x(0), v(0))/g(0, x(0), v(0)) - f_0(x(t), v(t))/g(0, x(0), v(0))
    static constexpr Real fOg{ 1 }; // f(0, x(0), v(0))/g(0, x(0), v(0)),
                                    // where g is the marker particle distribution

    explicit Particle() noexcept = default;
    explicit Particle(Vector const &vel, Real const pos_x) noexcept : vel{ vel }, pos_x{ pos_x } {}

    // pretty print
    //
    template <class CharT, class Traits>
    friend decltype(auto) operator<<(std::basic_ostream<CharT, Traits> &os, Particle const &ptl)
    {
        return os << '{' << ptl.vel << ", " << '{' << ptl.pos_x << '}' << ", " << '{' << ptl.f
                  << ", " << ptl.w << '}' << '}';
    }
};
PIC1D_END_NAMESPACE

#endif /* Particle_h */
