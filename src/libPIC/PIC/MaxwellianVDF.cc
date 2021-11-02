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
    vth1_eq       = std::sqrt(desc.beta1) * c * std::abs(desc.Oc) / desc.op;
    T2OT1_eq      = desc.T2_T1;
    sqrt_T2OT1_eq = std::sqrt(T2OT1_eq);
    vth1_eq_cubed = vth1_eq * vth1_eq * vth1_eq;
    //
    marker_vth1_eq       = vth1_eq * std::sqrt(desc.marker_temp_ratio);
    marker_vth1_eq_cubed = marker_vth1_eq * marker_vth1_eq * marker_vth1_eq;
    //
    N_extent.loc          = N(domain_extent.min());
    N_extent.len          = N(domain_extent.max()) - N_extent.loc;
    m_Nrefcell_div_Ntotal = (N(+0.5) - N(-0.5)) / N_extent.len;
}

auto MaxwellianVDF::eta(CurviCoord const &pos) const noexcept -> Real
{
    auto const cos = std::cos(geomtr.xi() * geomtr.D1() * pos.q1);
    return 1 / (T2OT1_eq + (1 - T2OT1_eq) * cos * cos);
}
auto MaxwellianVDF::T2OT1(CurviCoord const &pos) const noexcept -> Real
{
    return T2OT1_eq * eta(pos);
}
auto MaxwellianVDF::N(Real const q1) const noexcept -> Real
{
    if (geomtr.is_homogeneous())
        return q1;
    return std::atan(sqrt_T2OT1_eq * std::tan(geomtr.xi() * geomtr.D1() * q1)) / (sqrt_T2OT1_eq * geomtr.D1() * geomtr.xi());
}
auto MaxwellianVDF::q1(Real const N) const noexcept -> Real
{
    if (geomtr.is_homogeneous())
        return N;
    return std::atan(std::tan(sqrt_T2OT1_eq * geomtr.D1() * geomtr.xi() * N) / sqrt_T2OT1_eq) / (geomtr.xi() * geomtr.D1());
}

auto MaxwellianVDF::impl_n(CurviCoord const &pos) const -> Scalar
{
    constexpr Real n_eq = 1;
    return n_eq * eta(pos);
}
auto MaxwellianVDF::impl_nvv(CurviCoord const &pos) const -> Tensor
{
    Real const T1    = vth1_eq * vth1_eq;
    Real const T2OT1 = this->T2OT1(pos);
    Tensor     vv{ 1, T2OT1, T2OT1, 0, 0, 0 }; // field-aligned 2nd moment
    return geomtr.fac_to_cart(vv *= .5 * T1, pos) * Real{ impl_n(pos) };
}

auto MaxwellianVDF::f_common(Vector const &v, Real const T2OT1) noexcept -> Real
{
    // note that vel = {v1, v2, v3}/vth1
    // f0(x1, x2, x3) = exp(-x1^2)/√π * exp(-(x2^2 + x3^2)/(T2/T1))/(π T2/T1)
    //
    Real const f1 = std::exp(-v.x * v.x) * M_2_SQRTPI * .5;
    Real const x2 = v.y * v.y + v.z * v.z;
    Real const f2 = std::exp(-x2 / T2OT1) / (M_PI * T2OT1);
    return f1 * f2;
}
auto MaxwellianVDF::f0(Vector const &vel, CurviCoord const &pos) const noexcept -> Real
{
    return Real{ impl_n(pos) } * f_common(geomtr.cart_to_fac(vel, pos) / vth1_eq, T2OT1(pos)) / vth1_eq_cubed;
}
auto MaxwellianVDF::g0(Vector const &vel, CurviCoord const &pos) const noexcept -> Real
{
    return Real{ impl_n(pos) } * f_common(geomtr.cart_to_fac(vel, pos) / marker_vth1_eq, T2OT1(pos)) / marker_vth1_eq_cubed;
}

auto MaxwellianVDF::impl_emit(unsigned long const n) const -> std::vector<Particle>
{
    std::vector<Particle> ptls(n);
    for (auto &ptl : ptls)
        ptl = emit();
    return ptls;
}
auto MaxwellianVDF::impl_emit() const -> Particle
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
auto MaxwellianVDF::load() const -> Particle
{
    // position
    //
    CurviCoord const pos{ q1(bit_reversed<2>() * N_extent.len + N_extent.loc) };

    // velocity in field-aligned frame (Hu et al., 2010, doi:10.1029/2009JA015158)
    //
    Real const phi1 = bit_reversed<3>() * 2 * M_PI;                               // [0, 2pi]
    Real const v1   = std::sqrt(-std::log(uniform_real<100>())) * std::sin(phi1); // v_para
    //
    Real const phi2   = bit_reversed<5>() * 2 * M_PI; // [0, 2pi]
    Real const tmp_v2 = std::sqrt(-std::log(uniform_real<200>()) * T2OT1(pos));
    Real const v2     = std::cos(phi2) * tmp_v2; // in-plane v_perp
    Real const v3     = std::sin(phi2) * tmp_v2; // out-of-plane v_perp

    // velocity in Cartesian frame
    //
    Vector const vel = geomtr.fac_to_cart({ v1, v2, v3 }, pos) * marker_vth1_eq;

    return { vel, pos };
}
LIBPIC_END_NAMESPACE
