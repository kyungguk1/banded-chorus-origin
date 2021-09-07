/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "RelativisticMaxwellianVDF.h"
#include "RandomReal.h"
#include <cmath>

LIBPIC_BEGIN_NAMESPACE
RelativisticMaxwellianVDF::RelativisticMaxwellianVDF(BiMaxPlasmaDesc const &desc,
                                                     Geometry const        &geo,
                                                     Range const &domain_extent, Real c) noexcept
: RelativisticVDF{ geo, domain_extent, c }, desc{ desc }
{ // parameter check is assumed to be done already
    vth1       = std::sqrt(desc.beta1) * c * std::abs(desc.Oc) / desc.op;
    vth1_cubed = vth1 * vth1 * vth1;
    g2         = c * c / ((c - desc.Vd) * (c + desc.Vd));
    gd         = std::sqrt(g2);
    xd         = desc.Vd / vth1;
}

auto RelativisticMaxwellianVDF::n00c2(Real pos_x) const -> Scalar
{
    return n0(pos_x) * (c2 + .5 * vth1 * vth1 * (.5 + desc.T2_T1));
}
auto RelativisticMaxwellianVDF::P0Om0(Real pos_x) const -> Tensor
{
    Real const dP1    = std::sqrt(vth1 * vth1 / c2 * (0.75 + desc.T2_T1 / 2));
    Real const dP2    = std::sqrt(vth1 * vth1 / c2 * (0.25 + desc.T2_T1));
    Real const P1     = (1 - dP1) * (1 + dP1);
    Real const P2     = (1 - dP2) * (1 + dP2);
    Real const factor = n0(pos_x) * vth1 * vth1 / 2;
    return { factor * P1, factor * P2 * desc.T2_T1, factor * P2 * desc.T2_T1, 0, 0, 0 };
}

auto RelativisticMaxwellianVDF::impl_nU0(Real pos_x) const -> FourVector
{
    Real const ngc = n(pos_x) * gd * c;       // time component
    Real const ngV = n(pos_x) * desc.Vd * gd; // space components
    return geomtr.fac2cart({ ngc, { ngV, 0, 0 } });
}
auto RelativisticMaxwellianVDF::impl_nuv0(Real pos_x) const -> FourTensor
{
    Scalar const n00c2 = this->n00c2(pos_x);
    Tensor const P0Om0 = this->P0Om0(pos_x);
    Vector const Vd    = { desc.Vd, 0, 0 };
    Tensor const VV    = { Vd.x * Vd.x, 0, 0, 0, 0, 0 };

    // in field-aligned lab frame
    Scalar const ED = g2 * (n00c2 + P0Om0.xx * Vd.x / c2);
    Vector const MD = g2 / c2 * (*n00c2 + P0Om0.xx) * Vd;
    Tensor const uv = P0Om0 + g2 / c2 * (*n00c2 + P0Om0.xx) * VV;

    return { ED, geomtr.fac2cart(MD * c), geomtr.fac2cart(uv) };
}

inline auto RelativisticMaxwellianVDF::f0_comoving(const Vector &u0) const noexcept -> Real
{
    // note that u0 = {γv1, γv2, γv3}/vth1 in co-moving frame
    // f0(u1, u2, u3) = exp(-u1^2)/√π * exp(-(u2^2 + u3^2)/(T2/T1))/(π T2/T1)
    //
    Real const para = u0.x * u0.x;
    Real const f1   = std::exp(-para) * M_2_SQRTPI * .5;
    Real const perp = u0.y * u0.y + u0.z * u0.z;
    Real const f2   = std::exp(-perp / desc.T2_T1) / (M_PI * desc.T2_T1);
    return f1 * f2;
}
auto RelativisticMaxwellianVDF::f0_lab(Vector const &u) const noexcept -> Real
{
    // note that u = {γv1, γv2, γv3}/vth1 in lab frame
    // f(u1, u2, u3) = f0(γd(u1 - γu Vd), u2, u3)
    //
    Real const c2 = this->c2 / (vth1 * vth1);
    Real const gu = std::sqrt(1 + dot(u, u) / c2);
    Real const ux = gd * (u.x - gu * xd);
    return f0_comoving({ ux, u.y, u.z });
}

auto RelativisticMaxwellianVDF::impl_emit() const -> RelativisticParticle
{
    auto ptl = load();

    // rescale
    //
    ptl.g_vel *= vth1;
    ptl.pos_x *= domain_extent.len;
    ptl.pos_x += domain_extent.loc;

    // delta-f parameters
    //
    // ptl.fOg = ptl.f/ptl.g0(ptl);
    static_assert(Particle::PSD::fOg == 1.0, "f and g should be identical");
    ptl.psd = { f0(ptl), 1 };
    if (desc.scheme == ParticleScheme::delta_f) {
        ptl.psd.w = desc.initial_weight;
        ptl.psd.f /= 1 - ptl.psd.w / ptl.psd.fOg;
    }

    return ptl;
}
auto RelativisticMaxwellianVDF::load() const -> RelativisticParticle
{
    // position
    //
    Real const pos_x = bit_reversed<2>(); // [0, 1]

    // velocity in field-aligned co-moving frame (Hu et al., 2010, doi:10.1029/2009JA015158)
    //
    Real const phi1 = bit_reversed<3>() * 2 * M_PI;                               // [0, 2pi]
    Real const u1   = std::sqrt(-std::log(uniform_real<100>())) * std::sin(phi1); // γv_para
    //
    Real const phi2 = bit_reversed<5>() * 2 * M_PI; // [0, 2pi]
    Real const _u2  = std::sqrt(-std::log(uniform_real<200>()) * desc.T2_T1);
    Real const u2   = std::cos(phi2) * _u2; // in-plane γv_perp
    Real const u3   = std::sin(phi2) * _u2; // out-of-plane γv_perp

    // Lorentz transformation to lab frame
    //
    Real const   c2    = this->c2 / (vth1 * vth1);
    Real const   gu    = std::sqrt(1 + (u1 * u1 + u2 * u2 + u3 * u3) / c2);
    Vector const g_vel = { gd * (u1 + gu * xd), u2, u3 };

    return { geomtr.fac2cart(g_vel), pos_x, std::sqrt(1 + dot(g_vel, g_vel) / c2) };
}
LIBPIC_END_NAMESPACE
