/*
 * Copyright (c) 2021-2022, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC/RelativisticVDF.h>

LIBPIC_NAMESPACE_BEGIN(1)
/// Relativistic bi-Maxwellian velocity distribution function
/// \details
/// f(u1, u2) = n*exp(-x1^2 -x2^2)/(π^3/2 vth1^3 T2/T1),
///
/// where u = γv, x1 = u1/vth1, x2 = u2/(vth1*√(T2/T1))), and
/// T2 and T1 are temperatures in the directions perpendicular and
/// parallel to the background magnetic field direction, respectively.
///
class RelativisticMaxwellianVDF : public RelativisticVDF<RelativisticMaxwellianVDF> {
    using Super = RelativisticVDF<RelativisticMaxwellianVDF>;

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
    RelativisticMaxwellianVDF(BiMaxPlasmaDesc const &desc, Geometry const &geo, Range const &domain_extent, Real c);

    // VDF interfaces
    //
    [[nodiscard]] inline decltype(auto) impl_plasma_desc(Badge<Super>) const noexcept { return (this->desc); }

    [[nodiscard]] auto impl_n(Badge<Super>, CurviCoord const &pos) const -> Scalar;
    [[nodiscard]] auto impl_nV(Badge<Super>, CurviCoord const &) const -> CartVector
    {
        return { 0, 0, 0 };
    }
    [[nodiscard]] auto impl_nuv(Badge<Super>, CurviCoord const &pos) const -> FourCartTensor;

    [[nodiscard]] Real impl_Nrefcell_div_Ntotal(Badge<Super>) const { return m_Nrefcell_div_Ntotal; }
    [[nodiscard]] Real impl_f(Badge<Super>, Particle const &ptl) const { return f0(ptl); }

    [[nodiscard]] auto impl_emit(Badge<Super>, unsigned long) const -> std::vector<Particle>;
    [[nodiscard]] auto impl_emit(Badge<Super>) const -> Particle;

    // equilibrium physical distribution function
    //
    [[nodiscard]] Real f0(FourCartVector const &gcgvel, CurviCoord const &pos) const noexcept;
    [[nodiscard]] Real f0(Particle const &ptl) const noexcept { return f0(ptl.gcgvel, ptl.pos); }

    // marker particle distribution function
    //
    [[nodiscard]] Real g0(FourCartVector const &gcgvel, CurviCoord const &pos) const noexcept;
    [[nodiscard]] Real g0(Particle const &ptl) const noexcept { return g0(ptl.gcgvel, ptl.pos); }

private:
    [[nodiscard]] Real eta(CurviCoord const &pos) const noexcept;
    [[nodiscard]] Real T2OT1(CurviCoord const &pos) const noexcept;
    [[nodiscard]] Real N(Real q1) const noexcept;
    [[nodiscard]] Real q1(Real N) const noexcept;

    // total energy density; weakly relativistic limit
    [[nodiscard]] auto n00c2(CurviCoord const &pos) const -> Scalar;
    // pressure tensor; weakly relativistic limit
    [[nodiscard]] auto P0Om0(CurviCoord const &pos) const -> MFATensor;

    [[nodiscard]] Particle load() const;

    // velocity is normalized by vth1 and shifted to drifting plasma frame
    [[nodiscard]] static auto f_common(MFAVector const &g_vel, Real T2OT1, Real denom) noexcept -> Real;
};
LIBPIC_NAMESPACE_END(1)
