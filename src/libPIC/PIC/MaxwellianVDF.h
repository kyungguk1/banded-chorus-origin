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
/// f(v1, v2) = exp(-(x1 - xd)^2 -x2^2)/(π^3/2 vth1^3 T2/T1),
///
/// where x1 = v1/vth1, xd = vd/vth1, x2 = v2/(vth1*√(T2/T1))), and
/// T2 and T1 are temperatures in directions perpendicular and
/// parallel to the background magnetic field direction, respectively.
///
class MaxwellianVDF : public VDF<MaxwellianVDF> {
    friend VDF<MaxwellianVDF>;

    BiMaxPlasmaDesc desc;
    Real            vth1;  //!< Parallel thermal speed.
    Real            T2OT1; //!< Temperature anisotropy, T2/T1.
    Real            vth1_cubed;
    // marker psd parallel thermal speed
    Real marker_vth1;
    Real marker_vth1_cubed;

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

    [[nodiscard]] Scalar impl_n0(Real) const
    {
        constexpr Real n0 = 1;
        return n0;
    }
    [[nodiscard]] Vector impl_nV0(Real pos_x) const
    {
        return geomtr.fac2cart({ desc.Vd, 0, 0 }) * Real{ n0(pos_x) };
    }
    [[nodiscard]] Tensor impl_nvv0(Real pos_x) const
    {
        Real   xd = desc.Vd / vth1;
        Tensor vv{ 1 + 2 * xd * xd, T2OT1, T2OT1, 0, 0, 0 }; // field-aligned 2nd moment
        return geomtr.fac2cart(vv *= .5 * vth1 * vth1) * Real{ n0(pos_x) };
    }

    [[nodiscard]] Real impl_weight(Particle const &ptl) const { return (ptl.psd.real_f - f0(ptl)) / ptl.psd.marker; }

    [[nodiscard]] Particle impl_emit() const;
    [[nodiscard]] Particle load() const;

    // velocity is normalized by vth1
    [[nodiscard]] Real f0(Vector const &vel, Real xd) const noexcept;
    [[nodiscard]] Real f0(Vector const &vel) const noexcept { return f0(vel, desc.Vd / vth1); }
    [[nodiscard]] Real g0(Vector const &vel) const noexcept { return f0(vel, desc.Vd / marker_vth1); }

public:
    // equilibrium physical distribution function
    // f0(x1, x2, x3) = exp(-(x1 - xd)^2)/√π * exp(-(x2^2 + x3^2)/(T2/T1))/(π T2/T1)
    //
    [[nodiscard]] Real f0(Particle const &ptl) const noexcept
    {
        return f0(geomtr.cart2fac(ptl.vel) / vth1) * Real{ n0(ptl.pos_x) } / vth1_cubed;
    }

    // marker particle distribution function
    //
    [[nodiscard]] Real g0(Particle const &ptl) const noexcept
    {
        return g0(geomtr.cart2fac(ptl.vel) / marker_vth1) * Real{ n0(ptl.pos_x) } / marker_vth1_cubed;
    }
};
LIBPIC_END_NAMESPACE
