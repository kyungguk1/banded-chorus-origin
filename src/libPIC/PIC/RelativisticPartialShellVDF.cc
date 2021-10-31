/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "RelativisticPartialShellVDF.h"
#include "RandomReal.h"
#include <cmath>
#include <iterator>
#include <stdexcept>

LIBPIC_BEGIN_NAMESPACE
RelativisticPartialShellVDF::RelativisticPartialShellVDF(PartialShellPlasmaDesc const &desc, Geometry const &geo,
                                                         Range const &domain_extent, Real c) noexcept
: RelativisticVDF{ geo, domain_extent, c }, desc{ desc }
{
    vth             = std::sqrt(desc.beta) * c * std::abs(desc.Oc) / desc.op;
    auto const vth2 = vth * vth;
    vth_cubed       = vth2 * vth;
    {
        auto const xs  = desc.vs / vth;
        auto const xs2 = xs * xs;
        Ab             = .5 * (xs * std::exp(-xs2) + 2 / M_2_SQRTPI * (.5 + xs2) * std::erfc(-xs));
        Bz             = 2 / M_2_SQRTPI * std::tgamma(1 + .5 * desc.zeta) / std::tgamma(1.5 + .5 * desc.zeta);
        T_by_vth2      = .5 / Ab * (xs * (2.5 + xs2) * std::exp(-xs2) + 2 / M_2_SQRTPI * (0.75 + xs2 * (3 + xs2)) * std::erfc(-xs));
        TT_by_vth4     = .5 / Ab * (xs * (8.25 + xs2 * (7 + xs2)) * std::exp(-xs2) + 2 / M_2_SQRTPI * (1.875 + xs2 * (11.25 + 7.5 * xs2 + xs2 * xs2)) * std::erfc(-xs));
    }
    marker_vth       = vth * std::sqrt(desc.marker_temp_ratio);
    marker_vth_cubed = marker_vth * marker_vth * marker_vth;

    { // initialize q1 integral table
        N_extent.loc          = N_of_q1(domain_extent.min());
        N_extent.len          = N_of_q1(domain_extent.max()) - N_extent.loc;
        m_Nrefcell_div_Ntotal = (N_of_q1(+0.5) - N_of_q1(-0.5)) / N_extent.len;
        m_q1_of_N
            = init_integral_table(&RelativisticPartialShellVDF::N_of_q1, this, N_extent, domain_extent);
    }
    { // initialize velocity integral table
        constexpr Real t_max    = 5;
        auto const     xs       = desc.vs / marker_vth;
        Real const     t_min    = -(xs < t_max ? xs : t_max);
        Range const    x_extent = { t_min + xs, t_max - t_min };
        Fv_extent.loc           = Fv_of_x(x_extent.min());
        Fv_extent.len           = Fv_of_x(x_extent.max()) - Fv_extent.loc;
        m_x_of_Fv
            = init_integral_table(&RelativisticPartialShellVDF::Fv_of_x, this, Fv_extent, x_extent);
    }
    { // initialize pitch angle integral table
        constexpr auto accuracy_goal = 10;
        auto const     ph_max        = std::acos(std::pow(10, -accuracy_goal / Real(desc.zeta + 1)));
        auto const     ph_min        = -ph_max;
        Range const    a_extent      = { ph_min + M_PI_2, ph_max - ph_min };
        Fa_extent.loc                = Fa_of_a(a_extent.min());
        Fa_extent.len                = Fa_of_a(a_extent.max()) - Fa_extent.loc;
        m_a_of_Fa
            = init_integral_table(&RelativisticPartialShellVDF::Fa_of_a, this, Fa_extent, a_extent);
    }
}
auto RelativisticPartialShellVDF::init_integral_table(Real (RelativisticPartialShellVDF::*f_of_x)(Real) const noexcept,
                                                      RelativisticPartialShellVDF const *self, Range const f_extent, Range const x_extent)
    -> std::map<Real, Real>
{
    std::map<Real, Real> table;
    table.insert_or_assign(end(table), f_extent.min(), x_extent.min());
    constexpr long n_samples    = 5000;
    constexpr long n_subsamples = 100;
    auto const     df           = f_extent.len / n_samples;
    auto const     dx           = x_extent.len / (n_samples * n_subsamples);
    Real           x            = x_extent.min();
    Real           f_current    = (self->*f_of_x)(x);
    for (long i = 1; i < n_samples; ++i) {
        Real const f_target = Real(i) * df + f_extent.min();
        while (f_current < f_target)
            f_current = (self->*f_of_x)(x += dx);
        table.insert_or_assign(end(table), f_current, x);
    }
    table.insert_or_assign(end(table), f_extent.max(), x_extent.max());
    return table;
}

auto RelativisticPartialShellVDF::eta(CurviCoord const &pos) const noexcept -> Real
{
    if (desc.zeta == 0)
        return 1;
    auto const cos = std::cos(geomtr.xi() * geomtr.D1() * pos.q1);
    return std::pow(cos, desc.zeta);
}
auto RelativisticPartialShellVDF::int_cos_zeta(unsigned const zeta, Real const x) noexcept -> Real
{
    if (zeta == 0)
        return x;
    if (zeta == 1)
        return std::sin(x);
    return std::pow(std::cos(x), zeta - 1) * std::sin(x) / zeta + Real(zeta - 1) / zeta * int_cos_zeta(zeta - 2, x);
}
auto RelativisticPartialShellVDF::N_of_q1(Real const q1) const noexcept -> Real
{
    if (geomtr.is_homogeneous())
        return q1;
    return int_cos_zeta(desc.zeta, geomtr.xi() * geomtr.D1() * q1);
}
auto RelativisticPartialShellVDF::Fa_of_a(Real const alpha) const noexcept -> Real
{
    return int_cos_zeta(desc.zeta + 1, alpha - M_PI_2);
}
auto RelativisticPartialShellVDF::Fv_of_x(Real const v_by_vth) const noexcept -> Real
{
    auto const xs = desc.vs / marker_vth;
    auto const t  = v_by_vth - xs;
    return -(t + 2 * xs) * std::exp(-t * t) + 2 / M_2_SQRTPI * (.5 + xs * xs) * std::erf(t);
}
auto RelativisticPartialShellVDF::linear_interp(std::map<Real, Real> const &table, Real const x) -> Real
{
    auto const ub = table.upper_bound(x);
    if (ub == end(table) || ub == begin(table))
        throw std::out_of_range{ __PRETTY_FUNCTION__ };
    auto const &[x_min, y_min] = *std::prev(ub);
    auto const &[x_max, y_max] = *ub;
    return (y_min * (x_max - x) + y_max * (x - x_min)) / (x_max - x_min);
}
auto RelativisticPartialShellVDF::q1_of_N(Real const N) const -> Real
{
    return linear_interp(m_q1_of_N, N);
}
auto RelativisticPartialShellVDF::x_of_Fv(Real const Fv) const -> Real
{
    return linear_interp(m_x_of_Fv, Fv);
}
auto RelativisticPartialShellVDF::a_of_Fa(Real const Fa) const -> Real
{
    return linear_interp(m_a_of_Fa, Fa);
}

auto RelativisticPartialShellVDF::impl_n(CurviCoord const &pos) const -> Scalar
{
    constexpr Real n_eq = 1;
    return n_eq * eta(pos);
}

auto RelativisticPartialShellVDF::n00c2(CurviCoord const &pos) const -> Scalar
{
    return impl_n(pos) * (c2 + .5 * T_by_vth2 * vth * vth);
}
auto RelativisticPartialShellVDF::P0Om0(CurviCoord const &pos) const -> Tensor
{
    Real const P1 = *impl_n(pos) * vth * vth / Real(3 + desc.zeta) * (T_by_vth2 - .5 * TT_by_vth4 * vth * vth / c2);
    Real const P2 = P1 * (1 + .5 * desc.zeta);
    return { P1, P2, P2, 0, 0, 0 };
}

auto RelativisticPartialShellVDF::impl_nuv(CurviCoord const &pos) const -> FourTensor
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

auto RelativisticPartialShellVDF::f_common(Vector const &u_by_vth, unsigned const zeta, Real const us_by_vth, Real const Ab, Real const Bz) noexcept -> Real
{
    // note that x0 = {γv1, γv2, γv3}/vth in co-moving frame
    //
    // f(x1, x2, x3) = exp(-(x - xs)^2)*sin^ζ(α)/(2π θ^3 A(xs) B(ζ))
    //
    auto const x  = std::sqrt(dot(u_by_vth, u_by_vth));
    auto const t  = x - us_by_vth;
    Real const fv = std::exp(-t * t) / Ab;
    auto const u  = u_by_vth.x / x;
    Real const fa = (zeta == 0 ? 1 : std::pow((1 - u) * (1 + u), .5 * zeta)) / Bz;
    return .5 * fv * fa / M_PI;
}
auto RelativisticPartialShellVDF::f0(Vector const &vel, CurviCoord const &pos) const noexcept -> Real
{
    return Real{ impl_n(pos) } * f_common(geomtr.cart_to_fac(vel, pos) / vth, desc.zeta, desc.vs / vth, Ab, Bz) / vth_cubed;
}
auto RelativisticPartialShellVDF::g0(Vector const &vel, CurviCoord const &pos) const noexcept -> Real
{
    auto const xs        = desc.vs / marker_vth;
    auto const marker_Ab = .5 * (xs * std::exp(-xs * xs) + 2 / M_2_SQRTPI * (.5 + xs * xs) * std::erfc(-xs));
    return Real{ impl_n(pos) } * f_common(geomtr.cart_to_fac(vel, pos) / marker_vth, desc.zeta, xs, marker_Ab, Bz) / marker_vth_cubed;
}

auto RelativisticPartialShellVDF::impl_weight(Particle const &ptl) const -> Real
{
    switch (desc.scheme) {
        case ParticleScheme::full_f:
            return ptl.psd.real_f / ptl.psd.marker;
        case ParticleScheme::delta_f:
            return (ptl.psd.real_f - f0(ptl)) / ptl.psd.marker;
    }
}
auto RelativisticPartialShellVDF::impl_emit(unsigned long const n) const -> std::vector<Particle>
{
    std::vector<Particle> ptls(n);
    for (auto &ptl : ptls)
        ptl = emit();
    return ptls;
}
auto RelativisticPartialShellVDF::impl_emit() const -> Particle
{
    Particle ptl = load();

    switch (desc.scheme) {
        case ParticleScheme::full_f:
            ptl.psd        = { 1, f0(ptl), g0(ptl) };
            ptl.psd.weight = impl_weight(ptl);
            break;
        case ParticleScheme::delta_f:
            ptl.psd = { desc.initial_weight, f0(ptl), g0(ptl) };
            ptl.psd.real_f += ptl.psd.weight * ptl.psd.marker; // f = f_0 + w*g
            break;
    }

    return ptl;
}
auto RelativisticPartialShellVDF::load() const -> Particle
{
    // position
    //
    CurviCoord const pos{ q1_of_N(bit_reversed<2>() * N_extent.len + N_extent.loc) };

    // velocity in field-aligned frame (Hu et al., 2010, doi:10.1029/2009JA015158)
    //
    Real const ph     = bit_reversed<3>() * 2 * M_PI; // [0, 2pi]
    Real const alpha  = a_of_Fa(bit_reversed<5>() * Fa_extent.len + Fa_extent.loc);
    Real const u_vth  = x_of_Fv(uniform_real<100>() * Fv_extent.len + Fv_extent.loc);
    Real const u1     = std::cos(alpha) * u_vth;
    Real const tmp_u2 = std::sin(alpha) * u_vth;
    Real const u2     = std::cos(ph) * tmp_u2;
    Real const u3     = std::sin(ph) * tmp_u2;

    // velocity in Cartesian frame
    //
    Vector const g_vel = geomtr.fac_to_cart({ u1, u2, u3 }, pos) * marker_vth;
    auto const   gamma = std::sqrt(1 + dot(g_vel, g_vel) / c2);

    return { g_vel, pos, gamma };
}
LIBPIC_END_NAMESPACE
