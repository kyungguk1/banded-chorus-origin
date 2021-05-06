/*
 * Copyright (c) 2019-2021, Kyungguk Min
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

#ifndef Particle_h
#define Particle_h

#include "../Macros.h"
#include "../Predefined.h"
#include "./Vector.h"

#include <array>
#include <limits>
#include <ostream>

PIC1D_BEGIN_NAMESPACE
/// particle initialized by VDF
///
struct SimulationParticle {
    static constexpr Real quiet_nan = std::numeric_limits<Real>::quiet_NaN();

    Vector vel{ quiet_nan };   //!< 3-component velocity vector
    Real   pos_x{ quiet_nan }; //!< x-component of position

    // for delta-f
    //
    Real f{ quiet_nan }; // f(0, x(0), v(0))
    Real w{ quiet_nan }; // f(0, x(0), v(0))/g(0, x(0), v(0)) - f_0(x(t), v(t))/g(0, x(0), v(0))
    static constexpr Real fOg{ 1 }; // f(0, x(0), v(0))/g(0, x(0), v(0)),
                                    // where g is the marker particle distribution

    explicit constexpr SimulationParticle() noexcept = default;
    constexpr SimulationParticle(Vector const &vel, Real const pos_x) noexcept
    : vel{ vel }, pos_x{ pos_x }
    {
    }

private:
    // pretty print
    //
    template <class CharT, class Traits>
    friend decltype(auto) operator<<(std::basic_ostream<CharT, Traits> &os,
                                     SimulationParticle const &         ptl)
    {
        return os << '{' << ptl.vel << ", " << '{' << ptl.pos_x << '}' << ", " << '{' << ptl.f
                  << ", " << ptl.w << '}' << '}';
    }
};
static_assert(sizeof(SimulationParticle) == 6 * sizeof(Real));
static_assert(alignof(SimulationParticle) == alignof(Real));

/// single relativistic particle
///
struct RelativisticParticle : private SimulationParticle {
    using Super = SimulationParticle;

    using Super::f;
    using Super::pos_x;
    using Super::w;
    Real gamma{ quiet_nan }; //!< relativistic factor
    long padding{ -1 };

    explicit constexpr RelativisticParticle() noexcept = default;
    explicit constexpr RelativisticParticle(SimulationParticle const &ptl, Real gamma) noexcept
    : Super{ ptl }, gamma{ gamma }
    {
    }

    /// get SimulationParticle object
    [[nodiscard]] constexpr SimulationParticle const &base() const noexcept { return *this; }
    //[[nodiscard]] constexpr SimulationParticle &      base() noexcept { return *this; }

    /// gamma * velocity
    [[nodiscard]] constexpr Vector const &g_vel() const noexcept { return Super::vel; }
    [[nodiscard]] constexpr Vector &      g_vel() noexcept { return Super::vel; }

    /// usual velocity
    [[nodiscard]] constexpr Vector vel() const noexcept { return g_vel() / gamma; }

private:
    // pretty print
    //
    template <class CharT, class Traits>
    friend decltype(auto) operator<<(std::basic_ostream<CharT, Traits> &os,
                                     RelativisticParticle const &       ptl)
    {
        return os << '{' << ptl.base() << ", " << ptl.gamma << '}';
    }
};
static_assert(sizeof(RelativisticParticle) == 8 * sizeof(Real));
static_assert(alignof(RelativisticParticle) == alignof(Real));

/// single non-relativistic particle
///
struct NonrelativisticParticle : private SimulationParticle {
    using Super = SimulationParticle;

    using Super::f;
    using Super::pos_x;
    using Super::vel;
    using Super::w;
    std::array<long, 2> padding{ -1, -1 };

    explicit constexpr NonrelativisticParticle() noexcept = default;
    explicit constexpr NonrelativisticParticle(SimulationParticle const &ptl) noexcept
    : Super{ ptl }
    {
    }

    /// get SimulationParticle object
    [[nodiscard]] constexpr SimulationParticle const &base() const noexcept { return *this; }
    //[[nodiscard]] constexpr SimulationParticle &      base() noexcept { return *this; }

private:
    // pretty print
    //
    template <class CharT, class Traits>
    friend decltype(auto) operator<<(std::basic_ostream<CharT, Traits> &os,
                                     NonrelativisticParticle const &    ptl)
    {
        return os << ptl.base();
    }
};
static_assert(sizeof(NonrelativisticParticle) == 8 * sizeof(Real));
static_assert(alignof(NonrelativisticParticle) == alignof(Real));
PIC1D_END_NAMESPACE

#endif /* Particle_h */
