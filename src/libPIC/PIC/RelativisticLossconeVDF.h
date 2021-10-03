/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC/RelativisticVDF.h>
#include <map>

LIBPIC_BEGIN_NAMESPACE
/// Relativistic loss-cone velocity distribution function
/// \details
/// f(u1, u2) = n*exp(-x1^2)/(π^3/2 vth1 vth2^2) * (exp(-x2^2) - exp(-x2^2/β))/(1 - β),
///
/// where u = γv, x1 = u1/vth1, x2 = v2/vth2.
/// The effective temperature in the perpendicular direction is 2*T2/vth2^2 = 1 + β
///
class RelativisticLossconeVDF : public RelativisticVDF<RelativisticLossconeVDF> {
    friend RelativisticVDF<RelativisticLossconeVDF>;

public:
    struct RejectionSampler { // rejection sampler
        RejectionSampler() noexcept = default;
        explicit RejectionSampler(Real beta /*must not be 1*/);
        [[nodiscard]] inline Real sample() const noexcept;
        // ratio of the target to the proposed distributions
        [[nodiscard]] inline Real fOg(Real x) const noexcept;
        //
        static constexpr Real Delta{ 0 };     //!< Δ parameter.
        Real                  beta;           //!< β parameter.
        Real                  alpha;          //!< thermal spread of of the proposed distribution
        Real                  M;              //!< the ratio f(x_pk)/g(x_pk)
        static constexpr Real a_offset = 0.3; //!< optimal value for
                                              //!< thermal spread of the proposed distribution
    };

private:
    LossconePlasmaDesc desc;
    //
    Real beta_eq;         // loss-cone beta at the equator.
    Real vth1_eq;         //!< Parallel thermal speed at the equator.
    Real xth2_eq_squared; //!< The ratio of vth2^2 to vth1^2 at the equator.
    Real vth1_eq_cubed;
    // marker psd parallel thermal speed
    Real marker_vth1_eq;
    Real marker_vth1_eq_cubed;
    //
    Range                N_extent;
    Real                 m_Nrefcell_div_Ntotal;
    std::map<Real, Real> m_q1ofN;

public:
    /// Construct a loss-cone distribution
    /// \note Necessary parameter check is assumed to be done already.
    /// \param desc A BiMaxPlasmaDesc object.
    /// \param geo A geometry object.
    /// \param domain_extent Spatial domain extent.
    /// \param c Light speed. A positive real.
    ///
    RelativisticLossconeVDF(LossconePlasmaDesc const &desc, Geometry const &geo, Range const &domain_extent, Real c) noexcept;

private:
    [[nodiscard]] decltype(auto) impl_plasma_desc() const noexcept { return (this->desc); }

    [[nodiscard]] inline Real eta(CurviCoord const &pos) const noexcept;
    [[nodiscard]] inline Real eta_b(CurviCoord const &pos) const noexcept;
    [[nodiscard]] inline Real xth2_squared(CurviCoord const &pos) const noexcept;
    [[nodiscard]] inline Real beta(CurviCoord const &pos) const noexcept;
    [[nodiscard]] inline Real N(Real q1) const noexcept;
    [[nodiscard]] inline Real q1(Real N) const;

    // total energy density; weakly relativistic limit
    [[nodiscard]] inline Scalar n00c2(CurviCoord const &pos) const;
    // pressure tensor; weakly relativistic limit
    [[nodiscard]] inline Tensor P0Om0(CurviCoord const &pos) const;

    [[nodiscard]] Scalar     impl_n(CurviCoord const &pos) const;
    [[nodiscard]] Vector     impl_nV(CurviCoord const &) const { return { 0, 0, 0 }; }
    [[nodiscard]] FourTensor impl_nuv(CurviCoord const &pos) const;

    [[nodiscard]] Real impl_weight(Particle const &ptl) const { return (ptl.psd.real_f - f0(ptl)) / ptl.psd.marker; }

    [[nodiscard]] std::vector<Particle> impl_emit(unsigned long) const;
    [[nodiscard]] Particle              impl_emit() const;
    [[nodiscard]] Particle              load() const;

    // velocity is normalized by vth1
    [[nodiscard]] inline static Real f_common(Vector const &g_vel, Real xth2_squared, Real losscone_beta) noexcept;

public:
    // equilibrium physical distribution function
    //
    [[nodiscard]] Real f0(Vector const &g_vel, CurviCoord const &pos) const noexcept;
    [[nodiscard]] Real f0(Particle const &ptl) const noexcept { return f0(ptl.g_vel, ptl.pos); }

    // marker particle distribution function
    //
    [[nodiscard]] Real g0(Vector const &g_vel, CurviCoord const &pos) const noexcept;
    [[nodiscard]] Real g0(Particle const &ptl) const noexcept { return g0(ptl.g_vel, ptl.pos); }
};
LIBPIC_END_NAMESPACE
