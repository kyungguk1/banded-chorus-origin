/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC/VDF.h>
#include <map>

LIBPIC_BEGIN_NAMESPACE
/// Partial shell velocity distribution function
/// \details
/// f(v1, v2) = exp(-(x - xs)^2)*sin^ζ(α)/(2π θ^3 A(xs) B(ζ)),
///
/// where x = v/θ;
/// A(b) = (1/2) * (b exp(-b^2) + √π (1/2 + b^2) erfc(-b));
/// B(ζ) = √π Γ(1 + ζ/2)/Γ(1.5 + ζ/2).
///
class PartialShellVDF : public VDF<PartialShellVDF> {
    friend VDF<PartialShellVDF>;

    PartialShellPlasmaDesc desc;
    //
    Real vth;
    Real vth_cubed;
    Real Ab; // normalization factor for velocity distribution (for physical distribution)
    Real Bz; // normalization factor for pitch angle distribution
    Real T1; // parallel temperature, i.e., second velocity moment
    // marker psd parallel thermal speed
    Real marker_vth;
    Real marker_vth_cubed;
    //
    Range N_extent;
    Real  m_Nrefcell_div_Ntotal;
    Range Fv_extent;
    Range Fa_extent;
    //
    std::map<Real, Real> m_q1_of_N;
    std::map<Real, Real> m_x_of_Fv;
    std::map<Real, Real> m_a_of_Fa;

public:
    /// Construct a partial shell distribution
    /// \note Necessary parameter check is assumed to be done already.
    /// \param desc A PartialShellPlasmaDesc object.
    /// \param geo A geometry object.
    /// \param domain_extent Spatial domain extent.
    /// \param c Light speed. A positive real.
    ///
    PartialShellVDF(PartialShellPlasmaDesc const &desc, Geometry const &geo, Range const &domain_extent, Real c) noexcept;

private:
    [[nodiscard]] decltype(auto) impl_plasma_desc() const noexcept { return (this->desc); }

    [[nodiscard]] inline static Real          int_cos_zeta(unsigned zeta, Real x) noexcept;
    [[nodiscard]] static std::map<Real, Real> init_integral_table(Real (PartialShellVDF::*f_of_x)(Real) const noexcept,
                                                                  PartialShellVDF const *self, Range f_extent, Range x_extent);
    [[nodiscard]] static Real                 linear_interp(std::map<Real, Real> const &table, Real x);

    [[nodiscard]] inline Real eta(CurviCoord const &pos) const noexcept;
    [[nodiscard]] Real        N_of_q1(Real) const noexcept;
    [[nodiscard]] Real        Fa_of_a(Real) const noexcept;
    [[nodiscard]] Real        Fv_of_x(Real) const noexcept;
    [[nodiscard]] inline Real q1_of_N(Real) const;
    [[nodiscard]] inline Real x_of_Fv(Real) const;
    [[nodiscard]] inline Real a_of_Fa(Real) const;

    [[nodiscard]] Scalar impl_n(CurviCoord const &pos) const;
    [[nodiscard]] Vector impl_nV(CurviCoord const &) const { return { 0, 0, 0 }; }
    [[nodiscard]] Tensor impl_nvv(CurviCoord const &pos) const;

    [[nodiscard]] Real impl_weight(Particle const &ptl) const;

    [[nodiscard]] std::vector<Particle> impl_emit(unsigned long) const;
    [[nodiscard]] Particle              impl_emit() const;
    [[nodiscard]] Particle              load() const;

    // velocity is normalized by vth
    [[nodiscard]] inline static Real f_common(Vector const &vel_by_vth, unsigned zeta, Real vs_by_vth, Real Ab, Real Bz) noexcept;

public:
    // equilibrium physical distribution function
    //
    [[nodiscard]] Real f0(Vector const &vel, CurviCoord const &pos) const noexcept;
    [[nodiscard]] Real f0(Particle const &ptl) const noexcept { return f0(ptl.vel, ptl.pos); }

    // marker particle distribution function
    //
    [[nodiscard]] Real g0(Vector const &vel, CurviCoord const &pos) const noexcept;
    [[nodiscard]] Real g0(Particle const &ptl) const noexcept { return g0(ptl.vel, ptl.pos); }
};
LIBPIC_END_NAMESPACE
