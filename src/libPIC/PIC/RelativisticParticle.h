/*
 * Copyright (c) 2021, Kyungguk Min
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
/// single relativistic particle
///
struct RelativisticParticle {
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

    Vector g_vel{ quiet_nan }; //!< gamma * velocity, i.e., relativistic momentum
    Real   pos_x{ quiet_nan }; //!< x-component of position
    PSD    psd{};
    Real   gamma{ quiet_nan }; //!< relativistic factor; g = √(1 + g_vel^2/c^2)
    long   id{ -1 };           //!< particle identifier

    [[nodiscard]] Vector vel() const noexcept { return g_vel / gamma; } //!< Usual velocity

    RelativisticParticle() noexcept = default;
    RelativisticParticle(Vector const &g_vel, Real pos_x, Real gamma) noexcept
    : g_vel{ g_vel }, pos_x{ pos_x }, gamma{ gamma }, id{ next_id() }
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
    friend decltype(auto) operator<<(std::basic_ostream<CharT, Traits> &os, RelativisticParticle const &ptl)
    {
        return os << '{' << ptl.g_vel << ", " << '{' << ptl.pos_x << '}' << ", "
                  << '{' << ptl.psd.weight << ", " << ptl.psd.real_f << ", " << ptl.psd.marker << '}' << ", "
                  << ptl.gamma << ", " << ptl.id << '}';
    }
};
static_assert(sizeof(RelativisticParticle) == 9 * sizeof(RelativisticParticle::Real));
static_assert(alignof(RelativisticParticle) == alignof(Vector));
static_assert(std::is_standard_layout_v<RelativisticParticle>);
LIBPIC_END_NAMESPACE
