/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC/RelativisticVDF.h>

LIBPIC_BEGIN_NAMESPACE
/// Relativistic loss-cone velocity distribution function
/// \details
/// f0(u1, u2) = n0*exp(-x1^2)/(π^3/2 vth1 vth2^2) * ((1 - Δ*β)*exp(-x2^2) - (1 - / Δ)*exp(-x2^2/β))/(1 - β),
///
/// where u = γv, x1 = u1/vth1, x2 = v2/vth2.
/// The effective temperature in the perpendicular direction is 2*T2/vth2^2 = 1 + (1 - Δ)*β
///
/// In the lab frame,
/// f(u1, u2) = f0(γd(u1 - γVd), u2),
/// where γd = 1/√(1 - Vd^2/c^2) and γ = 1/√(1 - v^2/c^2).
///
class RelativisticLossconeVDF : public RelativisticVDF<RelativisticLossconeVDF> {
    friend RelativisticVDF<RelativisticLossconeVDF>;

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
    Real               vth1;       //!< Parallel thermal speed.
    Real               vth1_cubed; //!< vth1^3.
    Real               gd;         //!< γd.
    Real               g2;         //!< γd^2.
    Real               u2factor;   //!< 1 + (1 - Δ)β.
    Real               u4factor;   //!< 1 + (1 - Δ)(1 + β)β.
    Real               vth2;       //!< vth2
    // marker psd parallel thermal speed
    Real marker_vth1;
    Real marker_vth1_cubed;

public:
    /// Construct a loss-cone distribution
    /// \note Necessary parameter check is assumed to be done already.
    /// \param desc A BiMaxPlasmaDesc object.
    /// \param geo A geometry object.
    /// \param domain_extent Spatial domain extent.
    /// \param c Light speed. A positive real.
    ///
    RelativisticLossconeVDF(LossconePlasmaDesc const &desc, Geometry const &geo, Range const &domain_extent, Real c) noexcept;

    // number density in lab frame
    [[nodiscard]] static constexpr Real n_lab(Real) { return 1; }
    // number density in co-moving frame
    [[nodiscard]] Real n_comoving(Real pos_x) const { return n_lab(pos_x) / gd; }
    // total energy density in co-moving frame; non-relativistic limit
    [[nodiscard]] Scalar n00c2(Real pos_x) const;
    // pressure tensor in field-aligned co-moving frame; non-relativistic limit
    [[nodiscard]] Tensor P0Om0(Real pos_x) const;

private:
    [[nodiscard]] decltype(auto) impl_plasma_desc() const noexcept { return (this->desc); }

    [[nodiscard]] Scalar impl_n0(Real pos_x) const { return n_lab(pos_x); }
    [[nodiscard]] Vector impl_nV0(Real pos_x) const
    {
        return geomtr.fac2cart({ n_lab(pos_x) * desc.Vd, 0, 0 });
    }
    [[nodiscard]] FourTensor impl_nuv0(Real pos_x) const;

    [[nodiscard]] Real impl_weight(Particle const &ptl) const { return (ptl.psd.real_f - f0(ptl)) / ptl.psd.marker; }

    [[nodiscard]] Particle impl_emit() const;
    [[nodiscard]] Particle load() const;

    // velocity is normalized by vth1
    [[nodiscard]] Real f0_comoving(Vector const &g_vel) const noexcept;
    [[nodiscard]] Real f0_lab(Vector const &g_vel, Real denom) const noexcept;
    [[nodiscard]] Real f0_lab(Vector const &g_vel) const noexcept { return f0_lab(g_vel, vth1); }
    [[nodiscard]] Real g0_lab(Vector const &g_vel) const noexcept { return f0_lab(g_vel, marker_vth1); }

public:
    // equilibrium physical distribution function
    //
    [[nodiscard]] Real f0(Particle const &ptl) const noexcept
    {
        return f0_lab(geomtr.cart2fac(ptl.g_vel) / vth1) * n_comoving(ptl.pos_x) / vth1_cubed;
    }

    // marker particle distribution function
    //
    [[nodiscard]] Real g0(Particle const &ptl) const noexcept
    {
        return g0_lab(geomtr.cart2fac(ptl.g_vel) / marker_vth1) * n_comoving(ptl.pos_x) / marker_vth1_cubed;
    }
};
LIBPIC_END_NAMESPACE