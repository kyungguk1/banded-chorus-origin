/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "MaxwellianVDF.h"
#include "RandomReal.h"
#include <cmath>

LIBPIC_BEGIN_NAMESPACE
MaxwellianVDF::MaxwellianVDF(BiMaxPlasmaDesc const &desc, Geometry const &geo,
                             Range const &domain_extent, Real c) noexcept
: VDF{ geo, domain_extent }, desc{ desc }
{ // parameter check is assumed to be done already
    vth1       = std::sqrt(desc.beta1) * c * std::abs(desc.Oc) / desc.op;
    T2OT1      = desc.T2_T1;
    vth1_cubed = vth1 * vth1 * vth1;
    //
    marker_vth1       = vth1 * std::sqrt(desc.marker_temp_ratio);
    marker_vth1_cubed = marker_vth1 * marker_vth1 * marker_vth1;
}

auto MaxwellianVDF::f0(Vector const &v, Real const xd) const noexcept -> Real
{
    // note that vel = {v1, v2, v3}/vth1
    // f0(x1, x2, x3) = exp(-(x1 - xd)^2)/√π * exp(-(x2^2 + x3^2)/(T2/T1))/(π T2/T1)
    //
    Real const x1_xd      = v.x - xd;
    Real const f1         = std::exp(-x1_xd * x1_xd) * M_2_SQRTPI * .5;
    Real const x2_squared = v.y * v.y + v.z * v.z;
    Real const f2         = std::exp(-x2_squared / T2OT1) / (M_PI * T2OT1);
    return f1 * f2;
}

auto MaxwellianVDF::impl_emit() const -> Particle
{
    Particle ptl = load();

    // rescale
    //
    ptl.vel *= marker_vth1;
    ptl.pos_x *= domain_extent.len;
    ptl.pos_x += domain_extent.loc;

    // delta-f parameters
    //
    switch (desc.scheme) {
        case ParticleScheme::full_f:
            ptl.psd = { 1, -1, -1 };
            break;
        case ParticleScheme::delta_f:
            ptl.psd = { desc.initial_weight, f0(ptl), g0(ptl) };
            ptl.psd.real_f += ptl.psd.weight * ptl.psd.marker; // f = f_0 + w*g
            break;
    }

    return ptl;
}
auto MaxwellianVDF::load() const -> Particle
{
    // position
    //
    Real const pos_x = bit_reversed<2>(); // [0, 1]

    // velocity in field-aligned frame (Hu et al., 2010, doi:10.1029/2009JA015158)
    //
    Real const phi1 = bit_reversed<3>() * 2 * M_PI;                               // [0, 2pi]
    Real const v1   = std::sqrt(-std::log(uniform_real<100>())) * std::sin(phi1); // v_para
    //
    Real const phi2 = bit_reversed<5>() * 2 * M_PI; // [0, 2pi]
    Real const _v2  = std::sqrt(-std::log(uniform_real<200>()) * T2OT1);
    Real const v2   = std::cos(phi2) * _v2; // in-plane v_perp
    Real const v3   = std::sin(phi2) * _v2; // out-of-plane v_perp

    // velocity in Cartesian frame
    //
    Vector const vel = geomtr.fac2cart({ v1 + desc.Vd / marker_vth1, v2, v3 });

    return Particle{ vel, pos_x };
}
LIBPIC_END_NAMESPACE