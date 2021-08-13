/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <Geometry.h>
#include <PlasmaDesc.h>
#include <Predefined.h>
#include <Utility/Range.h>
#include <Utility/Scalar.h>
#include <Utility/Tensor.h>
#include <Utility/Vector.h>
#include <VDF/BitReversedPattern.h>
#include <VDF/Particle.h>
#include <common-config.h>

#include <memory>
#include <random>
#include <vector>

COMMON_BEGIN_NAMESPACE
/// Base class for velocity distribution function
///
class VDF {
protected:
    Geometry const geomtr;
    Range          domain_extent;

    VDF(Geometry const &geo, Range const &domain_extent)
    noexcept : geomtr{ geo }, domain_extent{ domain_extent }
    {
    }

public:
    virtual ~VDF() = default;

    /// Sample a single particle following the marker particle distribution, g0.
    [[nodiscard]] virtual Particle emit() const = 0;

    /// Sample N particles following the marker particle distribution, g0.
    /// \note This is a fallback; every subclass should implement this.
    [[nodiscard]] virtual std::vector<Particle> emit(unsigned n) const;

    /// Zero velocity moment at the given position, \<1\>_0(x).
    [[nodiscard]] virtual Scalar n0(Real pos_x) const = 0;

    /// First velocity moment at the given position, \<v\>_0(x).
    [[nodiscard]] virtual Vector nV0(Real pos_x) const = 0;

    /// Second velocity moment at the given position, \<vv\>_0(x).
    [[nodiscard]] virtual Tensor nvv0(Real) const = 0;

    /// Calculate the change of PSD
    /// \details Given a particle at some time instant, t,
    /// it calculate the change of PSD from the initial one:
    ///         1 - f_0(x(t), v(t))/f(0, x(0), v(0))
    ///
    [[nodiscard]] virtual Real delta_f(Particle const &ptl) const = 0;

    /// Calculate delta-f weighting factor
    [[nodiscard]] Real weight(Particle const &ptl) const
    {
        // f(0, x(0), v(0))/g(0, x(0), v(0)) - f_0(x(t), v(t))/g(0, x(0), v(0))
        // where g is the marker particle distribution
        //
        return ptl.psd.fOg * delta_f(ptl);
    }

private:
    // given a rng object, returns a real number, (0, 1), following a uniform distribution
    template <class URBG> [[nodiscard]] static Real uniform_real(URBG &g) noexcept
    {
        using uniform_t                   = std::uniform_real_distribution<>;
        static constexpr Real         eps = 1e-15;
        thread_local static uniform_t uniform{ eps, 1 - eps };
        return uniform(g);
    }

public:
    /// Returns a real number (0, 1) following a uniform distribution
    /// \tparam seed A seed for a random number generator.
    ///
    template <unsigned seed> [[nodiscard]] static Real uniform_real() noexcept
    { // seed must be passed as a template parameter
        static_assert(seed > 0, "seed has to be a positive number");
        thread_local static std::mt19937 g{ seed };
        return uniform_real(g);
    }
    /// Returns a real number (0, 1) following a uniform distribution
    /// \tparam base Base prime number for BitReversedPattern.
    ///
    template <unsigned base> [[nodiscard]] static Real bit_reversed() noexcept
    {
        static_assert(base > 0, "base has to be a positive number");
        thread_local static BitReversedPattern<base> g{};
        return uniform_real(g);
    }
};
COMMON_END_NAMESPACE
