/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC/RelativisticVDF.h>

LIBPIC_BEGIN_NAMESPACE
/// Relativistic bi-Maxwellian velocity distribution function
/// \details
/// f(u1, u2) = n*exp(-x1^2 -x2^2)/(π^3/2 vth1^3 T2/T1),
///
/// where u = γv, x1 = u1/vth1, x2 = u2/(vth1*√(T2/T1))), and
/// T2 and T1 are temperatures in the directions perpendicular and
/// parallel to the background magnetic field direction, respectively.
///
class RelativisticMaxwellianVDF : public RelativisticVDF<RelativisticMaxwellianVDF> {
    friend RelativisticVDF<RelativisticMaxwellianVDF>;

    BiMaxPlasmaDesc desc;
    //
    Real vth1_eq; //!< Parallel thermal speed.
    Real vth1_eq_cubed;
    Real T2OT1_eq; //!< Temperature anisotropy, T2/T1, at the equator.
    Real sqrt_T2OT1_eq;
    // marker psd parallel thermal speed
    Real marker_vth1_eq;
    Real marker_vth1_eq_cubed;
    //
    Range N_extent;
    Real  m_Nrefcell_div_Ntotal;

public:
    /// Construct a relativistic bi-Maxwellian distribution
    /// \note Necessary parameter check is assumed to be done already.
    /// \param desc A BiMaxPlasmaDesc object.
    /// \param geo A geometry object.
    /// \param domain_extent Spatial domain extent.
    /// \param c Light speed. A positive real.
    ///
    RelativisticMaxwellianVDF(BiMaxPlasmaDesc const &desc, Geometry const &geo, Range const &domain_extent, Real c) noexcept;

private:
    [[nodiscard]] decltype(auto) impl_plasma_desc() const noexcept { return (this->desc); }

    [[nodiscard]] inline Real eta(CurviCoord const &pos) const noexcept;
    [[nodiscard]] inline Real T2OT1(CurviCoord const &pos) const noexcept;
    [[nodiscard]] inline Real N(Real q1) const noexcept;
    [[nodiscard]] inline Real q1(Real N) const noexcept;

    // total energy density; weakly relativistic limit
    [[nodiscard]] inline Scalar n00c2(CurviCoord const &pos) const;
    // pressure tensor; weakly relativistic limit
    [[nodiscard]] inline Tensor P0Om0(CurviCoord const &pos) const;

    [[nodiscard]] Scalar     impl_n(CurviCoord const &pos) const;
    [[nodiscard]] Vector     impl_nV(CurviCoord const &) const { return { 0, 0, 0 }; }
    [[nodiscard]] FourTensor impl_nuv(CurviCoord const &pos) const;

    [[nodiscard]] Real impl_weight(Particle const &ptl) const;

    [[nodiscard]] std::vector<Particle> impl_emit(unsigned long) const;
    [[nodiscard]] Particle              impl_emit() const;
    [[nodiscard]] Particle              load() const;

    // velocity is normalized by vth1
    [[nodiscard]] inline static Real f_common(Vector const &g_vel, Real T2OT1) noexcept;

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
