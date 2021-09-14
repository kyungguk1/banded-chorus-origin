/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "RelativisticLossconeVDF.h"
#include "RandomReal.h"
#include <cmath>
#include <stdexcept>

LIBPIC_BEGIN_NAMESPACE
RelativisticLossconeVDF::RelativisticLossconeVDF(LossconePlasmaDesc const &desc, Geometry const &geo,
                                                 Range const &domain_extent, Real c) noexcept
: RelativisticVDF{ geo, domain_extent, c }, desc{ desc }
{ // parameter check is assumed to be done already
    Real const Delta = desc.Delta;
    Real const beta  = [beta = desc.beta]() noexcept { // avoid beta == 1
        constexpr Real eps  = 1e-5;
        Real const     diff = beta - 1;
        return beta + (std::abs(diff) < eps ? std::copysign(eps, diff) : 0);
    }();
    rs         = RejectionSampler{ Delta, beta };
    vth1       = std::sqrt(desc.beta1) * c * std::abs(desc.Oc) / desc.op;
    vth1_cubed = vth1 * vth1 * vth1;
    g2         = c * c / ((c - desc.Vd) * (c + desc.Vd));
    gd         = std::sqrt(g2);
    xd         = desc.Vd / vth1;
    //
    u2factor = 1 + (1 - Delta) * beta;
    u4factor = 1 + (1 - Delta) * (1 + beta) * beta;
    vth2     = vth1 * std::sqrt(desc.T2_T1 / u2factor);
}

auto RelativisticLossconeVDF::n00c2(Real pos_x) const -> Scalar
{
    return n_comoving(pos_x) * (c2 + .5 * vth1 * vth1 * (.5 + desc.T2_T1));
}
auto RelativisticLossconeVDF::P0Om0(Real pos_x) const -> Tensor
{
    Real const dP1    = std::sqrt(vth1 * vth1 / c2 * (0.75 + desc.T2_T1 / 2));
    Real const dP2    = std::sqrt(vth1 * vth1 / c2 * (0.25 + desc.T2_T1 * u4factor / (u2factor * u2factor)));
    Real const P1     = (1 - dP1) * (1 + dP1);
    Real const P2     = (1 - dP2) * (1 + dP2);
    Real const factor = n_comoving(pos_x) * vth1 * vth1 / 2;
    return { factor * P1, factor * P2 * desc.T2_T1, factor * P2 * desc.T2_T1, 0, 0, 0 };
}

auto RelativisticLossconeVDF::impl_nuv0(Real pos_x) const -> FourTensor
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

inline auto RelativisticLossconeVDF::f0_comoving(const Vector &u0) const noexcept -> Real
{
    // note that x0 = {γv1, γv2, γv3}/vth1 in co-moving frame
    //
    //                                  ((1 - Δ*β)*exp(-(x2^2 + x3^2)/xth2^2) - (1 - Δ)*exp(-(x2^2 + x3^2)/(β*xth2^2)))
    // f0(x1, x2, x3) = exp(-x1^2)/√π * -------------------------------------------------------------------------------
    //                                                             (π * xth2^2 * (1 - β))
    //
    Real const para = u0.x * u0.x;
    Real const f1   = std::exp(-para) * M_2_SQRTPI * .5;
    Real const xth2 = vth2 / vth1;
    Real const f2   = [D = rs.Delta, b = rs.beta,
                     t2 = (u0.y * u0.y + u0.z * u0.z) / (xth2 * xth2),
                     de = M_PI * (xth2 * xth2) * (1 - rs.beta)]() noexcept {
        return ((1 - D * b) * std::exp(-t2) - (1 - D) * std::exp(-t2 / b)) / de;
    }();
    return f1 * f2;
}
auto RelativisticLossconeVDF::f0_lab(Vector const &u) const noexcept -> Real
{
    // note that u = {γv1, γv2, γv3}/vth1 in lab frame
    // f(u1, u2, u3) = f0(γd(u1 - γu Vd), u2, u3)
    //
    Real const c2 = this->c2 / (vth1 * vth1);
    Real const gu = std::sqrt(1 + dot(u, u) / c2);
    Real const ux = gd * (u.x - gu * xd);
    return f0_comoving({ ux, u.y, u.z });
}

auto RelativisticLossconeVDF::impl_emit() const -> Particle
{
    Particle ptl = load();

    // rescale
    //
    ptl.g_vel *= vth1;
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
auto RelativisticLossconeVDF::load() const -> Particle
{
    // position
    //
    Real const pos_x = bit_reversed<2>(); // [0, 1]

    // velocity in field-aligned co-moving frame (Hu et al., 2010, doi:10.1029/2009JA015158)
    //
    Real const phi1 = bit_reversed<3>() * 2 * M_PI;                               // [0, 2pi]
    Real const u1   = std::sqrt(-std::log(uniform_real<100>())) * std::sin(phi1); // v_para
    //
    Real const phi2 = bit_reversed<5>() * 2 * M_PI; // [0, 2pi]
    Real const _u2  = rs.sample() * vth2 / vth1;
    Real const u2   = std::cos(phi2) * _u2; // in-plane γv_perp
    Real const u3   = std::sin(phi2) * _u2; // out-of-plane γv_perp

    // Lorentz transformation to lab frame
    //
    Real const   c2    = this->c2 / (vth1 * vth1);
    Real const   gu    = std::sqrt(1 + (u1 * u1 + u2 * u2 + u3 * u3) / c2);
    Vector const g_vel = { gd * (u1 + gu * xd), u2, u3 };

    return { geomtr.fac2cart(g_vel), pos_x, std::sqrt(1 + dot(g_vel, g_vel) / c2) };
}

// MARK: - RejectionSampler
//
RelativisticLossconeVDF::RejectionSampler::RejectionSampler(Real const Delta,
                                                            Real const beta /*must not be 1*/)
: Delta{ Delta }, beta{ beta }
{
    constexpr Real eps = 1e-5;
    if (std::abs(1 - Delta) < eps) { // Δ == 1
        alpha = 1;
        M     = 1;
    } else { // Δ != 1
        alpha               = (beta < 1 ? 1 : beta) + a_offset;
        auto const eval_xpk = [D = Delta, b = beta, a = alpha] {
            Real const det
                = -b / (1 - b) * std::log(((a - 1) * (1 - D * b) * b) / ((a - b) * (1 - D)));
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
