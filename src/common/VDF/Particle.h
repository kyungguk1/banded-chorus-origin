/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

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
    struct PSD {
        Real f{ quiet_nan }; // f(0, x(0), v(0))
        Real w{ quiet_nan }; // f(0, x(0), v(0))/g(0, x(0), v(0)) - f_0(x(t), v(t))/g(0, x(0), v(0))
        static constexpr Real fOg{ 1 }; // f(0, x(0), v(0))/g(0, x(0), v(0)),
                                        // where g is the marker particle distribution
    };

    Vector vel{ quiet_nan };   //!< 3-component velocity vector
    Real   pos_x{ quiet_nan }; //!< x-component of position
    PSD    psd{};
    long   id{ -1 }; //!< particle identifier

    Particle() noexcept = default;
    Particle(Vector const &vel, Real pos_x) noexcept : vel{ vel }, pos_x{ pos_x }, id{ next_id() }
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
        return os << '{' << ptl.vel << ", " << '{' << ptl.pos_x << '}' << ", " << '{' << ptl.psd.f
                  << ", " << ptl.psd.w << '}' << '}';
    }
};
static_assert(sizeof(Particle) == 8 * sizeof(Particle::Real));
static_assert(alignof(Particle) == alignof(Vector));
static_assert(std::is_standard_layout_v<Particle>);

/// single relativistic particle
///
struct RelativisticParticle {
    using Real = Particle::Real;
    using PSD  = Particle::PSD;

    Vector g_vel{ Particle::quiet_nan }; //!< gamma * velocity, i.e., relativistic momentum
    Real   pos_x{ Particle::quiet_nan }; //!< x-component of position
    PSD    psd{};
    Real   gamma{ Particle::quiet_nan }; //!< relativistic factor; g = √(1 + g_vel^2/c^2)
    [[nodiscard]] Vector vel() const noexcept { return g_vel / gamma; } //!< Usual velocity

    RelativisticParticle() noexcept = default;
    RelativisticParticle(Vector const &g_vel, Real pos_x, Real gamma) noexcept
    : g_vel{ g_vel }, pos_x{ pos_x }, gamma{ gamma }
    {
    }
    RelativisticParticle(Particle const &ptl, Real gamma) noexcept
    : RelativisticParticle{ ptl.vel, ptl.pos_x, gamma }
    {
    }

    // pretty print
    //
    template <class CharT, class Traits>
    friend decltype(auto) operator<<(std::basic_ostream<CharT, Traits> &os,
                                     RelativisticParticle const        &ptl)
    {
        return os << '{' << ptl.g_vel << ", " << '{' << ptl.pos_x << '}' << ", " << '{' << ptl.psd.f
                  << ", " << ptl.psd.w << '}' << ", " << ptl.gamma << '}';
    }
};
static_assert(sizeof(RelativisticParticle) == sizeof(Particle));
static_assert(alignof(RelativisticParticle) == alignof(Particle));
static_assert(std::is_standard_layout_v<RelativisticParticle>);
COMMON_END_NAMESPACE
