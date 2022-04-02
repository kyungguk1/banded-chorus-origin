/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC/VDF.h>

LIBPIC_BEGIN_NAMESPACE
/// Bi-Maxwellian velocity distribution function
/// \details
/// f(v1, v2) = exp(-x1^2 -x2^2)/(π^3/2 vth1^3 T2/T1),
///
/// where x1 = v1/vth1, x2 = v2/(vth1*√(T2/T1))), and
/// T2 and T1 are temperatures in directions perpendicular and
/// parallel to the background magnetic field direction, respectively.
///
class MaxwellianVDF : public VDF<MaxwellianVDF> {
    friend VDF<MaxwellianVDF>;

    BiMaxPlasmaDesc desc;
    //
    Real vth1_eq;  //!< Parallel thermal speed at the equator.
    Real T2OT1_eq; //!< Temperature anisotropy, T2/T1, at the equator.
    Real sqrt_T2OT1_eq;
    Real vth1_eq_cubed;
    // marker psd parallel thermal speed
    Real marker_vth1_eq;
    Real marker_vth1_eq_cubed;
    //
    Range N_extent;
    Real  m_Nrefcell_div_Ntotal;

public:
    /// Construct a bi-Maxwellian distribution
    /// \note Necessary parameter check is assumed to be done already.
    /// \param desc A BiMaxPlasmaDesc object.
    /// \param geo A geometry object.
    /// \param domain_extent Spatial domain extent.
    /// \param c Light speed. A positive real.
    ///
    MaxwellianVDF(BiMaxPlasmaDesc const &desc, Geometry const &geo, Range const &domain_extent, Real c) noexcept;

private:
    [[nodiscard]] decltype(auto) impl_plasma_desc() const noexcept { return (this->desc); }

    [[nodiscard]] inline Real eta(CurviCoord const &pos) const noexcept;
    [[nodiscard]] inline Real T2OT1(CurviCoord const &pos) const noexcept;
    [[nodiscard]] inline Real N(Real q1) const noexcept;
    [[nodiscard]] inline Real q1(Real N) const noexcept;

    [[nodiscard]] Scalar impl_n(CurviCoord const &pos) const;
    [[nodiscard]] Vector impl_nV(CurviCoord const &) const { return { 0, 0, 0 }; }
    [[nodiscard]] Tensor impl_nvv(CurviCoord const &pos) const;

    [[nodiscard]] Real impl_f0(Particle const &ptl) const { return f0(ptl); }

    [[nodiscard]] std::vector<Particle> impl_emit(unsigned long) const;
    [[nodiscard]] Particle              impl_emit() const;
    [[nodiscard]] Particle              load() const;

    // velocity is normalized by vth1
    [[nodiscard]] inline static Real f_common(Vector const &vel, Real T2OT1) noexcept;

public:
    // equilibrium physical distribution function
    // f0(x1, x2, x3) = exp(-x1^2)/√π * exp(-(x2^2 + x3^2)/(T2/T1))/(π T2/T1)
    //
    [[nodiscard]] Real f0(Vector const &vel, CurviCoord const &pos) const noexcept;
    [[nodiscard]] Real f0(Particle const &ptl) const noexcept { return f0(ptl.vel, ptl.pos); }

    // marker particle distribution function
    //
    [[nodiscard]] Real g0(Vector const &vel, CurviCoord const &pos) const noexcept;
    [[nodiscard]] Real g0(Particle const &ptl) const noexcept { return g0(ptl.vel, ptl.pos); }
};
LIBPIC_END_NAMESPACE
