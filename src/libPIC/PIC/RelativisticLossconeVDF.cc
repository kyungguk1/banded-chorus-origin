/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "RelativisticLossconeVDF.h"
#include "RandomReal.h"
#include <cmath>
#include <iterator>
#include <stdexcept>

LIBPIC_BEGIN_NAMESPACE
RelativisticLossconeVDF::RelativisticLossconeVDF(LossconePlasmaDesc const &desc, Geometry const &geo,
                                                 Range const &domain_extent, Real c) noexcept
: RelativisticVDF{ geo, domain_extent, c }, desc{ desc }
{
    beta_eq = [](Real const beta) noexcept { // avoid beta == 1 && beta == 0
        constexpr Real eps = 1e-5;
        if (beta < eps)
            return eps;
        if (Real const diff = beta - 1; std::abs(diff) < eps)
            return beta + std::copysign(eps, diff);
        return beta;
    }(desc.beta);
    //
    vth1_eq         = std::sqrt(desc.beta1) * c * std::abs(desc.Oc) / desc.op;
    xth2_eq_squared = desc.T2_T1 / (1 + beta_eq);
    vth1_eq_cubed   = vth1_eq * vth1_eq * vth1_eq;
    //
    marker_vth1_eq       = vth1_eq * std::sqrt(desc.marker_temp_ratio);
    marker_vth1_eq_cubed = marker_vth1_eq * marker_vth1_eq * marker_vth1_eq;
    //
    N_extent.loc          = N(domain_extent.min());
    N_extent.len          = N(domain_extent.max()) - N_extent.loc;
    m_Nrefcell_div_Ntotal = (N(+0.5) - N(-0.5)) / N_extent.len;
    //
    m_q1ofN.insert_or_assign(end(m_q1ofN), N_extent.min(), domain_extent.min());
    constexpr long n_samples    = 50000;
    constexpr long n_subsamples = 100;
    auto const     dN           = N_extent.len / n_samples;
    auto const     dq           = domain_extent.len / (n_samples * n_subsamples);
    Real           q1           = domain_extent.min();
    Real           N_current    = N(q1);
    for (long i = 1; i < n_samples; ++i) {
        Real const N_target = Real(i) * dN + N_extent.min();
        while (N_current < N_target)
            N_current = N(q1 += dq);
        m_q1ofN.insert_or_assign(end(m_q1ofN), N_current, q1);
    }
    m_q1ofN.insert_or_assign(end(m_q1ofN), N_extent.max(), domain_extent.max());
}

auto RelativisticLossconeVDF::eta(CurviCoord const &pos) const noexcept -> Real
{
    auto const cos = std::cos(geomtr.xi() * geomtr.D1() * pos.q1);
    return 1 / (xth2_eq_squared + (1 - xth2_eq_squared) * cos * cos);
}
auto RelativisticLossconeVDF::eta_b(CurviCoord const &pos) const noexcept -> Real
{
    auto const cos = std::cos(geomtr.xi() * geomtr.D1() * pos.q1);
    auto const tmp = beta_eq * xth2_eq_squared;
    return 1 / (tmp + (1 - tmp) * cos * cos);
}
auto RelativisticLossconeVDF::xth2_squared(CurviCoord const &pos) const noexcept -> Real
{
    return xth2_eq_squared * eta(pos);
}
auto RelativisticLossconeVDF::beta(CurviCoord const &pos) const noexcept -> Real
{
    auto const beta = beta_eq * eta_b(pos) / eta(pos);
    // avoid beta == 1
    constexpr Real eps = 1e-5;
    if (Real const diff = beta - 1; std::abs(diff) < eps)
        return beta + std::copysign(eps, diff);
    return beta;
}
auto RelativisticLossconeVDF::N(Real const q1) const noexcept -> Real
{
    if (geomtr.is_homogeneous())
        return q1;

    auto const xth2_eq      = std::sqrt(xth2_eq_squared);
    auto const sqrt_beta_eq = std::sqrt(beta_eq);
    auto const tan          = std::tan(geomtr.xi() * geomtr.D1() * q1);
    auto const tmp1         = std::atan(xth2_eq * tan) / (xth2_eq * geomtr.D1() * geomtr.xi());
    auto const tmp2         = std::atan(sqrt_beta_eq * xth2_eq * tan) / (sqrt_beta_eq * xth2_eq * geomtr.D1() * geomtr.xi());
    return (tmp1 - beta_eq * tmp2) / (1 - beta_eq);
}
auto RelativisticLossconeVDF::q1(Real const N) const -> Real
{
    auto const ub = m_q1ofN.upper_bound(N);
    if (ub == end(m_q1ofN) || ub == begin(m_q1ofN))
        throw std::out_of_range{ __PRETTY_FUNCTION__ };
    auto const lb = std::prev(ub);
    // linear interpolation
    return (lb->second * (ub->first - N) + ub->second * (N - lb->first)) / (ub->first - lb->first);
}

auto RelativisticLossconeVDF::impl_n(CurviCoord const &pos) const -> Scalar
{
    constexpr Real n_eq = 1;
    return n_eq * (eta(pos) - beta_eq * eta_b(pos)) / (1 - beta_eq);
}

auto RelativisticLossconeVDF::n00c2(CurviCoord const &pos) const -> Scalar
{
    Real const T2OT1 = (1 + beta(pos)) * xth2_squared(pos);
    return impl_n(pos) * (c2 + .5 * vth1_eq * vth1_eq * (.5 + T2OT1));
}
auto RelativisticLossconeVDF::P0Om0(CurviCoord const &pos) const -> Tensor
{
    Real const vth1_squared = vth1_eq * vth1_eq;
    Real const beta         = this->beta(pos);
    auto const u2factor     = 1 + beta;
    auto const u4factor     = 1 + (1 + beta) * beta;
    Real const T2OT1        = (1 + beta) * xth2_squared(pos);
    Real const dP1          = std::sqrt(vth1_squared / c2 * (0.75 + T2OT1 / 2));
    Real const dP2          = std::sqrt(vth1_squared / c2 * (0.25 + T2OT1 * u4factor / (u2factor * u2factor)));
    Real const P1           = (1 - dP1) * (1 + dP1);
    Real const P2           = (1 - dP2) * (1 + dP2);
    Real const factor       = *impl_n(pos) * vth1_squared / 2;
    return { factor * P1, factor * P2 * T2OT1, factor * P2 * T2OT1, 0, 0, 0 };
}

auto RelativisticLossconeVDF::impl_nuv(CurviCoord const &pos) const -> FourTensor
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

auto RelativisticLossconeVDF::f_common(Vector const &u, Real const xth2_squared, Real const beta) noexcept -> Real
{
    // note that x0 = {γv1, γv2, γv3}/vth1 in co-moving frame
    //
    //                                  (exp(-(x2^2 + x3^2)/xth2^2) - exp(-(x2^2 + x3^2)/(β*xth2^2)))
    // f0(x1, x2, x3) = exp(-x1^2)/√π * -------------------------------------------------------------
    //                                                     (π * xth2^2 * (1 - β))
    //
    Real const f1 = std::exp(-u.x * u.x) * M_2_SQRTPI * .5;
    Real const f2 = [D = 0, b = beta, t2 = (u.y * u.y + u.z * u.z) / xth2_squared,
                     de = M_PI * xth2_squared * (1 - beta)]() noexcept {
        return ((1 - D * b) * std::exp(-t2) - (1 - D) * std::exp(-t2 / b)) / de;
    }();
    return f1 * f2;
}
auto RelativisticLossconeVDF::f0(Vector const &vel, CurviCoord const &pos) const noexcept -> Real
{
    return Real{ impl_n(pos) } * f_common(geomtr.cart_to_fac(vel, pos) / vth1_eq, xth2_squared(pos), beta(pos)) / vth1_eq_cubed;
}
auto RelativisticLossconeVDF::g0(Vector const &vel, CurviCoord const &pos) const noexcept -> Real
{
    return Real{ impl_n(pos) } * f_common(geomtr.cart_to_fac(vel, pos) / marker_vth1_eq, xth2_squared(pos), beta(pos)) / marker_vth1_eq_cubed;
}

auto RelativisticLossconeVDF::impl_emit(unsigned long const n) const -> std::vector<Particle>
{
    std::vector<Particle> ptls(n);
    for (auto &ptl : ptls)
        ptl = emit();
    return ptls;
}
auto RelativisticLossconeVDF::impl_emit() const -> Particle
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
auto RelativisticLossconeVDF::load() const -> Particle
{
    // position
    //
    CurviCoord const pos{ q1(bit_reversed<2>() * N_extent.len + N_extent.loc) };

    // velocity in field-aligned co-moving frame (Hu et al., 2010, doi:10.1029/2009JA015158)
    //
    Real const phi1 = bit_reversed<3>() * 2 * M_PI;                               // [0, 2pi]
    Real const u1   = std::sqrt(-std::log(uniform_real<100>())) * std::sin(phi1); // v_para
    //
    Real const phi2   = bit_reversed<5>() * 2 * M_PI; // [0, 2pi]
    Real const tmp_u2 = RejectionSampler{ beta(pos) }.sample() * std::sqrt(xth2_squared(pos));
    Real const u2     = std::cos(phi2) * tmp_u2; // in-plane v_perp
    Real const u3     = std::sin(phi2) * tmp_u2; // out-of-plane v_perp

    // velocity in Cartesian frame
    //
    Vector const g_vel = geomtr.fac_to_cart({ u1, u2, u3 }, pos) * marker_vth1_eq;
    auto const   gamma = std::sqrt(1 + dot(g_vel, g_vel) / c2);

    return { g_vel, pos, gamma };
}

// MARK: - RejectionSampler
//
RelativisticLossconeVDF::RejectionSampler::RejectionSampler(Real const beta /*must not be 1*/)
: beta{ beta }
{
    constexpr Real eps = 1e-5;
    if (std::abs(1 - Delta) < eps) { // Δ == 1
        alpha = 1;
        M     = 1;
    } else { // Δ != 1
        alpha               = (beta < 1 ? 1 : beta) + a_offset;
        auto const eval_xpk = [D = Delta, b = beta, a = alpha] {
            Real const det = -b / (1 - b) * std::log(((a - 1) * (1 - D * b) * b) / ((a - b) * (1 - D)));
            return std::isfinite(det) && det > 0 ? std::sqrt(det) : 0;
        };
        Real const xpk = std::abs(1 - Delta * beta) < eps ? 0 : eval_xpk();
        M              = fOg(xpk);
    }
    if (!std::isfinite(M))
        throw std::runtime_error{ __PRETTY_FUNCTION__ };
}
auto RelativisticLossconeVDF::RejectionSampler::fOg(const Real x) const noexcept -> Real
{
    using std::exp;
    Real const x2 = x * x;
    Real const f  = ((1 - Delta * beta) * exp(-x2) - (1 - Delta) * exp(-x2 / beta)) / (1 - beta);
    Real const g  = exp(-x2 / alpha) / alpha;
    return f / g; // ratio of the target distribution to proposed distribution
}
auto RelativisticLossconeVDF::RejectionSampler::sample() const noexcept -> Real
{
    auto const vote = [this](Real const proposal) noexcept {
        Real const jury = uniform_real<300>() * M;
        return jury <= fOg(proposal);
    };
    auto const proposed = [a = this->alpha]() noexcept {
        return std::sqrt(-std::log(uniform_real<200>()) * a);
    };
    //
    Real sample;
    while (!vote(sample = proposed())) {}
    return sample;
}
LIBPIC_END_NAMESPACE
