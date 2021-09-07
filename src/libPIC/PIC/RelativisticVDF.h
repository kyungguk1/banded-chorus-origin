/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC/Config.h>
#include <PIC/FourTensor.h>
#include <PIC/FourVector.h>
#include <PIC/Geometry.h>
#include <PIC/Particle.h>
#include <PIC/PlasmaDesc.h>
#include <PIC/Predefined.h>
#include <PIC/Range.h>
#include <PIC/Scalar.h>

#include <vector>

LIBPIC_BEGIN_NAMESPACE
/// Base class for velocity distribution function
///
template <class Concrete> class RelativisticVDF {
    using Self = Concrete;

    [[nodiscard]] constexpr auto self() const noexcept { return static_cast<Self const *>(this); }
    [[nodiscard]] constexpr auto self() noexcept { return static_cast<Self *>(this); }

protected:
    Geometry geomtr;
    Range    domain_extent;
    Real     c;
    Real     c2;

    RelativisticVDF(Geometry const &geo, Range const &domain_extent, Real c) noexcept
    : geomtr{ geo }, domain_extent{ domain_extent }, c{ c }, c2{ c * c }
    {
    }

public:
    /// Plasma description associated with *this
    ///
    [[nodiscard]] decltype(auto) plasma_desc() const noexcept { return self()->impl_plasma_desc(); }

    /// Sample a single particle following the marker particle distribution, g0.
    /// \note Concrete subclass should provide impl_emit with the same signature.
    ///
    [[nodiscard]] RelativisticParticle emit() const { return self()->impl_emit(); }

    /// Sample N particles following the marker particle distribution, g0.
    [[nodiscard]] auto emit(unsigned n) const
    {
        std::vector<RelativisticParticle> ptls(n);
        for (auto &ptl : ptls)
            ptl = emit();
        return ptls;
    }

    /// Particle density-flux four vector, γ*{\<c\>_0(x), \<v\>_0(x)}.
    /// \details Here γ = 1/√(1 - Vd^2/c^2).
    /// \note Concrete subclass should provide impl_nU0 with the same signature.
    ///
    [[nodiscard]] FourVector nU0(Real pos_x) const { return self()->impl_nU0(pos_x); }

    /// Stress-energy four-tensor, \<ui*vj\>_0(x).
    /// \details Here u = {γc, γv} and v = {c, v}.
    /// \note Concrete subclass should provide impl_nuv0 with the same signature.
    ///
    [[nodiscard]] FourTensor nuv0(Real pos_x) const { return self()->impl_nuv0(pos_x); }

    /// Calculate the change of PSD
    /// \details Given a particle at some time instant, t,
    /// it calculate the change of PSD from the initial one:
    ///         1 - f_0(x(t), u(t))/f(0, x(0), u(0))
    ///
    /// \note Concrete subclass should provide impl_delta_f with the same signature.
    ///
    [[nodiscard]] Real delta_f(RelativisticParticle const &ptl) const
    {
        return self()->impl_delta_f(ptl);
    }

    /// Calculate delta-f weighting factor
    /// \details The weight is given by
    ///
    /// f(0, x(0), u(0))/g(0, x(0), u(0)) - f_0(x(t), u(t))/g(0, x(0), u(0))
    ///
    /// where g is the marker particle distribution.
    ///
    [[nodiscard]] Real weight(RelativisticParticle const &ptl) const
    {
        return ptl.psd.fOg * delta_f(ptl);
    }
};
LIBPIC_END_NAMESPACE
