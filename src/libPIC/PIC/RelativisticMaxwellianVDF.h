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
/// In the co-moving frame,
/// f0(u1, u2) = n0*exp(-x1^2 -x2^2)/(π^3/2 vth1^3 T2/T1),
///
/// where u = γv, x1 = u1/vth1, x2 = u2/(vth1*√(T2/T1))), and
/// T2 and T1 are temperatures in the directions perpendicular and
/// parallel to the background magnetic field direction, respectively.
///
/// In the lab frame,
/// f(u1, u2) = f0(γd(u1 - γVd), u2),
/// where γd = 1/√(1 - Vd^2/c^2) and γ = 1/√(1 - v^2/c^2).
///
class RelativisticMaxwellianVDF : public RelativisticVDF<RelativisticMaxwellianVDF> {
    friend RelativisticVDF<RelativisticMaxwellianVDF>;

    BiMaxPlasmaDesc desc;
    Real            vth1; //!< Parallel thermal speed.
    Real            vth1_cubed;
    Real            gd; //!< γd.
    Real            g2; //!< γd^2.
    Real            xd; //!< Vd/vth1

public:
    /// Construct a relativistic bi-Maxwellian distribution
    /// \note Necessary parameter check is assumed to be done already.
    /// \param desc A BiMaxPlasmaDesc object.
    /// \param geo A geometry object.
    /// \param domain_extent Spatial domain extent.
    /// \param c Light speed. A positive real.
    ///
    RelativisticMaxwellianVDF(BiMaxPlasmaDesc const &desc, Geometry const &geo, Range const &domain_extent, Real c) noexcept;

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
    [[nodiscard]] Real f0_lab(Vector const &g_vel) const noexcept;
    [[nodiscard]] Real g0_lab(Vector const &g_vel) const noexcept { return f0_lab(g_vel); }

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
        return g0_lab(geomtr.cart2fac(ptl.g_vel) / vth1) * n_comoving(ptl.pos_x) / vth1_cubed;
    }
};
LIBPIC_END_NAMESPACE
