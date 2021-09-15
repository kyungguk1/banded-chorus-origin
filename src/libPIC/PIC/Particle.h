/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC/Config.h>
#include <PIC/Vector.h>

#include <limits>
#include <ostream>
#include <type_traits>

LIBPIC_BEGIN_NAMESPACE
/// single particle description
///
struct Particle {
    using Real                      = double;
    static constexpr auto quiet_nan = std::numeric_limits<Real>::quiet_NaN();

    // for delta-f
    struct PSD {
        Real weight{ quiet_nan }; // f(0, x(0), v(0))/g(0, x(0), v(0)) - f_0(x(t), v(t))/g(0, x(0), v(0))
        Real real_f{ quiet_nan }; // f(0, x(0), v(0)), physical particle PSD
        Real marker{ quiet_nan }; // g(0, x(0), v(0)), simulation particle PSD

        constexpr PSD() noexcept = default;
        constexpr PSD(Real w, Real f, Real g) noexcept
        : weight{ w }, real_f{ f }, marker{ g } {}
    };

    Vector vel{ quiet_nan };   //!< 3-component velocity vector
    Real   pos_x{ quiet_nan }; //!< x-component of position
    PSD    psd{};
    long   id{ -1 }; //!< particle identifier

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
                  << '{' << ptl.psd.weight << ", " << ptl.psd.real_f << ", " << ptl.psd.marker << '}' << ", "
                  << ptl.id << '}';
    }
};
static_assert(sizeof(Particle) == 8 * sizeof(Particle::Real));
static_assert(alignof(Particle) == alignof(Vector));
static_assert(std::is_standard_layout_v<Particle>);
LIBPIC_END_NAMESPACE
