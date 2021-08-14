/*
 * Copyright (c) 2020-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC/VDF.h>

LIBPIC_BEGIN_NAMESPACE
/// Loss-cone velocity distribution function
/// \details
/// f(v1, v2) = exp(-(x1 - xd)^2)/(π^3/2 vth1 vth2^2) * ((1 - Δ*β)*exp(-x2^2) -
/// (1 - / Δ)*exp(-x2^2/β))/(1 - β),
///
/// where x1 = v1/vth1, xd = vd/vth1, x2 = v2/vth2.
/// The effective temperature in the perpendicular direction is 2*T2/vth2^2 = 1 + (1 - Δ)*β
///
class LossconeVDF : public VDF<LossconeVDF> {
    friend VDF<LossconeVDF>;

public:
    struct RejectionSampler { // rejection sampler
        RejectionSampler() noexcept = default;
        RejectionSampler(Real Delta, Real beta /*must not be 1*/);
        [[nodiscard]] inline Real sample() const noexcept;
        // ratio of the target to the proposed distributions
        [[nodiscard]] inline Real fOg(Real x) const noexcept;
        //
        Real Delta; //!< Δ parameter.
        Real beta;  //!< β parameter.
    public:
        Real                  alpha;          //!< thermal spread of of the proposed distribution
        Real                  M;              //!< the ratio f(x_pk)/g(x_pk)
        static constexpr Real a_offset = 0.3; //!< optimal value for
                                              //!< thermal spread of the proposed distribution
    };

private:
    LossconePlasmaDesc desc;
    RejectionSampler   rs;
    Real               vth1;         //!< Parallel thermal speed.
    Real               T2OT1;        //!< Temperature anisotropy, T2/T1.
    Real               xd;           //!< Parallel drift speed normalized to vth1.
    Real               xth2_squared; //!< The ratio of vth2^2 to vth1^2.
    Real               vth1_cubed;   //!< vth1^3.

public:
    /// Construct a loss-cone distribution
    /// \note Necessary parameter check is assumed to be done already.
    /// \param desc A BiMaxPlasmaDesc object.
    /// \param geo A geometry object.
    /// \param domain_extent Spatial domain extent.
    /// \param c Light speed. A positive real.
    ///
    LossconeVDF(LossconePlasmaDesc const &desc, Geometry const &geo, Range const &domain_extent,
                Real c) noexcept;

private:
    [[nodiscard]] Scalar impl_n0(Real) const
    {
        constexpr Real n0 = 1;
        return n0;
    }
    [[nodiscard]] Vector impl_nV0(Real pos_x) const
    {
        return geomtr.fac2cart({ xd * vth1, 0, 0 }) * Real{ n0(pos_x) };
    }
    [[nodiscard]] Tensor impl_nvv0(Real pos_x) const
    {
        Tensor vv{ 1 + 2 * xd * xd, T2OT1, T2OT1, 0, 0, 0 }; // field-aligned 2nd moment
        return geomtr.fac2cart(vv *= .5 * vth1 * vth1) * Real{ n0(pos_x) };
    }

    [[nodiscard]] Real impl_delta_f(Particle const &ptl) const { return 1 - f0(ptl) / ptl.psd.f; }

    [[nodiscard]] Particle impl_emit() const;
    [[nodiscard]] Particle load() const;

    // equilibrium physical distribution function
    //
    [[nodiscard]] Real f0(Vector const &vel) const noexcept;
    [[nodiscard]] Real f0(Particle const &ptl) const noexcept
    {
        return f0(geomtr.cart2fac(ptl.vel) / vth1) / vth1_cubed;
    }

    // marker particle distribution function
    //
    [[nodiscard]] Real g0(Vector const &vel) const noexcept { return f0(vel); }
    [[nodiscard]] Real g0(Particle const &ptl) const noexcept
    {
        return g0(geomtr.cart2fac(ptl.vel) / vth1) / vth1_cubed;
    }
};
LIBPIC_END_NAMESPACE
