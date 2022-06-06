/*
 * Copyright (c) 2021-2022, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "MaxwellianVDF.h"
#include "../RandomReal.h"
#include <algorithm>
#include <cmath>

LIBPIC_NAMESPACE_BEGIN(1)
RelativisticMaxwellianVDF::RelativisticMaxwellianVDF(BiMaxPlasmaDesc const &desc, Geometry const &geo, Range const &domain_extent, Real c)
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

auto RelativisticMaxwellianVDF::impl_n(Badge<Super>, CurviCoord const &pos) const -> Scalar
{
    constexpr Real n_eq = 1;
    return n_eq * eta(pos);
}

auto RelativisticMaxwellianVDF::n00c2(CurviCoord const &pos) const -> Scalar
{
    return this->n0(pos) * (c2 + .5 * vth1_eq * vth1_eq * (.5 + T2OT1(pos)));
}
auto RelativisticMaxwellianVDF::P0Om0(CurviCoord const &pos) const -> MFATensor
{
    Real const vth1_squared = vth1_eq * vth1_eq;
    Real const T2OT1        = this->T2OT1(pos);
    Real const dP1          = std::sqrt(vth1_squared / c2 * (0.75 + T2OT1 / 2));
    Real const dP2          = std::sqrt(vth1_squared / c2 * (0.25 + T2OT1));
    Real const P1           = (1 - dP1) * (1 + dP1);
    Real const P2           = (1 - dP2) * (1 + dP2);
    Real const factor       = *this->n0(pos) * vth1_squared / 2;
    return { factor * P1, factor * P2 * T2OT1, factor * P2 * T2OT1, 0, 0, 0 };
}
auto RelativisticMaxwellianVDF::impl_nuv(Badge<Super>, CurviCoord const &pos) const -> FourCartTensor
{
    Scalar const    n00c2 = this->n00c2(pos);
    MFATensor const P0Om0 = this->P0Om0(pos);
    MFAVector const Vd    = { 0, 0, 0 };
    MFATensor const VV    = { Vd.x * Vd.x, 0, 0, 0, 0, 0 };

    // in field-aligned lab frame
    constexpr auto  g2 = 1;
    Scalar const    ED = g2 * (n00c2 + P0Om0.xx * Vd.x / c2);
    MFAVector const MD = g2 / c2 * (*n00c2 + P0Om0.xx) * Vd;
    MFATensor const uv = P0Om0 + g2 / c2 * (*n00c2 + P0Om0.xx) * VV;

    return { ED, geomtr.mfa_to_cart(MD * c, pos), geomtr.mfa_to_cart(uv, pos) };
}

auto RelativisticMaxwellianVDF::f_common(MFAVector const &u0, Real const T2OT1, Real const denom) noexcept -> Real
{
    // note that the u0 is in co-moving frame and normalized by vth1
    // f0(u1, u2, u3) = exp(-u1^2)/√π * exp(-(u2^2 + u3^2)/(T2/T1))/(π T2/T1)
    //
    Real const f1   = std::exp(-u0.x * u0.x) * M_2_SQRTPI * .5;
    Real const perp = u0.y * u0.y + u0.z * u0.z;
    Real const f2   = std::exp(-perp / T2OT1) / (M_PI * T2OT1);
    return (f1 * f2) / denom;
}
auto RelativisticMaxwellianVDF::f0(FourCartVector const &gcgvel, CurviCoord const &pos) const noexcept -> Real
{
    // note that u = γ{v1, v2, v3} in lab frame, where γ = c/√(c^2 - v^2)
    auto const gcgv_mfa = geomtr.cart_to_mfa(gcgvel, pos);
    return Real{ this->n0(pos) } * f_common(gcgv_mfa.s / vth1_eq, T2OT1(pos), vth1_eq_cubed);
}
auto RelativisticMaxwellianVDF::g0(FourCartVector const &gcgvel, CurviCoord const &pos) const noexcept -> Real
{
    // note that u = γ{v1, v2, v3} in lab frame, where γ = c/√(c^2 - v^2)
    auto const gcgv_mfa = geomtr.cart_to_mfa(gcgvel, pos);
    return Real{ this->n0(pos) } * f_common(gcgv_mfa.s / marker_vth1_eq, T2OT1(pos), marker_vth1_eq_cubed);
}

auto RelativisticMaxwellianVDF::impl_emit(Badge<Super>, unsigned long const n) const -> std::vector<Particle>
{
    std::vector<Particle> ptls(n);
    std::generate(begin(ptls), end(ptls), [this] {
        return this->emit();
    });
    return ptls;
}
auto RelativisticMaxwellianVDF::impl_emit(Badge<Super>) const -> Particle
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
    Real const phi2 = bit_reversed<5>() * 2 * M_PI; // [0, 2pi]
    Real const tmp  = std::sqrt(-std::log(uniform_real<200>()) * T2OT1(pos));
    Real const u2   = std::cos(phi2) * tmp; // in-plane γv_perp
    Real const u3   = std::sin(phi2) * tmp; // out-of-plane γv_perp

    // boost from particle reference frame to co-moving frame
    auto const gcgv_mfa = lorentz_boost<-1>(FourMFAVector{ c, {} }, MFAVector{ u1, u2, u3 } * (marker_vth1_eq / c));

    return { geomtr.mfa_to_cart(gcgv_mfa, pos), pos };
}
LIBPIC_NAMESPACE_END(1)
