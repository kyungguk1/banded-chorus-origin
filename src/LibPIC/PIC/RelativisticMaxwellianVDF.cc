/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "RelativisticMaxwellianVDF.h"
#include "RandomReal.h"
#include <cmath>

LIBPIC_BEGIN_NAMESPACE
RelativisticMaxwellianVDF::RelativisticMaxwellianVDF(BiMaxPlasmaDesc const &desc, Geometry const &geo,
                                                     Range const &domain_extent, Real c) noexcept
: RelativisticVDF{ geo, domain_extent, c }, desc{ desc }
{ // parameter check is assumed to be done already
    vth1_eq       = std::sqrt(desc.beta1) * c * std::abs(desc.Oc) / desc.op;
    vth1_eq_cubed = vth1_eq * vth1_eq * vth1_eq;
    T2OT1_eq      = desc.T2_T1;
    sqrt_T2OT1_eq = std::sqrt(T2OT1_eq);
    //
    marker_vth1_eq       = vth1_eq * std::sqrt(desc.marker_temp_ratio);
    marker_vth1_eq_cubed = marker_vth1_eq * marker_vth1_eq * marker_vth1_eq;
    //
    N_extent.loc          = N(domain_extent.min());
    N_extent.len          = N(domain_extent.max()) - N_extent.loc;
    m_Nrefcell_div_Ntotal = (N(+0.5) - N(-0.5)) / N_extent.len;
}

auto RelativisticMaxwellianVDF::eta(CurviCoord const &pos) const noexcept -> Real
{
    auto const cos = std::cos(geomtr.xi() * geomtr.D1() * pos.q1);
    return 1 / (T2OT1_eq + (1 - T2OT1_eq) * cos * cos);
}
auto RelativisticMaxwellianVDF::T2OT1(CurviCoord const &pos) const noexcept -> Real
{
    return T2OT1_eq * eta(pos);
}
auto RelativisticMaxwellianVDF::N(Real const q1) const noexcept -> Real
{
    if (geomtr.is_homogeneous())
        return q1;
    return std::atan(sqrt_T2OT1_eq * std::tan(geomtr.xi() * geomtr.D1() * q1)) / (sqrt_T2OT1_eq * geomtr.D1() * geomtr.xi());
}
auto RelativisticMaxwellianVDF::q1(Real const N) const noexcept -> Real
{
    if (geomtr.is_homogeneous())
        return N;
    return std::atan(std::tan(sqrt_T2OT1_eq * geomtr.D1() * geomtr.xi() * N) / sqrt_T2OT1_eq) / (geomtr.xi() * geomtr.D1());
}

auto RelativisticMaxwellianVDF::impl_n(CurviCoord const &pos) const -> Scalar
{
    constexpr Real n_eq = 1;
    return n_eq * eta(pos);
}

auto RelativisticMaxwellianVDF::n00c2(CurviCoord const &pos) const -> Scalar
{
    return impl_n(pos) * (c2 + .5 * vth1_eq * vth1_eq * (.5 + T2OT1(pos)));
}
auto RelativisticMaxwellianVDF::P0Om0(CurviCoord const &pos) const -> Tensor
{
    Real const vth1_squared = vth1_eq * vth1_eq;
    Real const T2OT1        = this->T2OT1(pos);
    Real const dP1          = std::sqrt(vth1_squared / c2 * (0.75 + T2OT1 / 2));
    Real const dP2          = std::sqrt(vth1_squared / c2 * (0.25 + T2OT1));
    Real const P1           = (1 - dP1) * (1 + dP1);
    Real const P2           = (1 - dP2) * (1 + dP2);
    Real const factor       = *impl_n(pos) * vth1_squared / 2;
    return { factor * P1, factor * P2 * T2OT1, factor * P2 * T2OT1, 0, 0, 0 };
}
auto RelativisticMaxwellianVDF::impl_nuv(CurviCoord const &pos) const -> FourTensor
{
    Scalar const n00c2 = this->n00c2(pos);
    Tensor const P0Om0 = this->P0Om0(pos);
    Vector const Vd    = { 0, 0, 0 };
    Tensor const VV    = { Vd.x * Vd.x, 0, 0, 0, 0, 0 };

    // in field-aligned lab frame
    constexpr auto g2 = 1;
    Scalar const   ED = g2 * (n00c2 + P0Om0.xx * Vd.x / c2);
    Vector const   MD = g2 / c2 * (*n00c2 + P0Om0.xx) * Vd;
    Tensor const   uv = P0Om0 + g2 / c2 * (*n00c2 + P0Om0.xx) * VV;

    return { ED, geomtr.fac_to_cart(MD * c, pos), geomtr.fac_to_cart(uv, pos) };
}

auto RelativisticMaxwellianVDF::f_common(Vector const &u, Real const T2OT1) noexcept -> Real
{
    // note that u0 = {γv1, γv2, γv3}/vth1 in co-moving frame
    // f0(u1, u2, u3) = exp(-u1^2)/√π * exp(-(u2^2 + u3^2)/(T2/T1))/(π T2/T1)
    //
    Real const f1 = std::exp(-u.x * u.x) * M_2_SQRTPI * .5;
    Real const x2 = u.y * u.y + u.z * u.z;
    Real const f2 = std::exp(-x2 / T2OT1) / (M_PI * T2OT1);
    return f1 * f2;
}
auto RelativisticMaxwellianVDF::f0(Vector const &vel, CurviCoord const &pos) const noexcept -> Real
{
    return Real{ impl_n(pos) } * f_common(geomtr.cart_to_fac(vel, pos) / vth1_eq, T2OT1(pos)) / vth1_eq_cubed;
}
auto RelativisticMaxwellianVDF::g0(Vector const &vel, CurviCoord const &pos) const noexcept -> Real
{
    return Real{ impl_n(pos) } * f_common(geomtr.cart_to_fac(vel, pos) / marker_vth1_eq, T2OT1(pos)) / marker_vth1_eq_cubed;
}

auto RelativisticMaxwellianVDF::impl_emit(unsigned long const n) const -> std::vector<Particle>
{
    std::vector<Particle> ptls(n);
    for (auto &ptl : ptls)
        ptl = emit();
    return ptls;
}
auto RelativisticMaxwellianVDF::impl_emit() const -> Particle
{
    Particle ptl = load();

    switch (desc.scheme) {
        case ParticleScheme::full_f:
            ptl.psd        = { 1, f0(ptl), g0(ptl) };
            ptl.psd.weight = ptl.psd.real_f / ptl.psd.marker;
            break;
        case ParticleScheme::delta_f:
            ptl.psd = { desc.initial_weight, f0(ptl), g0(ptl) };
            ptl.psd.real_f += ptl.psd.weight * ptl.psd.marker; // f = f_0 + w*g
            break;
    }

    return ptl;
}
auto RelativisticMaxwellianVDF::load() const -> Particle
{
    // position
    //
    CurviCoord const pos{ q1(bit_reversed<2>() * N_extent.len + N_extent.loc) };

    // velocity in field-aligned co-moving frame (Hu et al., 2010, doi:10.1029/2009JA015158)
    //
    Real const phi1 = bit_reversed<3>() * 2 * M_PI;                               // [0, 2pi]
    Real const u1   = std::sqrt(-std::log(uniform_real<100>())) * std::sin(phi1); // γv_para
    //
    Real const phi2   = bit_reversed<5>() * 2 * M_PI; // [0, 2pi]
    Real const tmp_u2 = std::sqrt(-std::log(uniform_real<200>()) * T2OT1(pos));
    Real const u2     = std::cos(phi2) * tmp_u2; // in-plane γv_perp
    Real const u3     = std::sin(phi2) * tmp_u2; // out-of-plane γv_perp

    // velocity in Cartesian frame
    //
    Vector const g_vel = geomtr.fac_to_cart({ u1, u2, u3 }, pos) * marker_vth1_eq;
    auto const   gamma = std::sqrt(1 + dot(g_vel, g_vel) / c2);

    return { g_vel, pos, gamma };
}
LIBPIC_END_NAMESPACE
