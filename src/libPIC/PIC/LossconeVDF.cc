/*
 * Copyright (c) 2020-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "LossconeVDF.h"
#include "RandomReal.h"
#include <cmath>
#include <stdexcept>
#include <utility>

LIBPIC_BEGIN_NAMESPACE
LossconeVDF::LossconeVDF(Geometry const &geo, Range const &domain_extent,
                         LossconePlasmaDesc const &desc, Real c) noexcept
: VDF{ geo, domain_extent }, desc{ desc }
{ // parameter check is assumed to be done already
    Real const Delta = desc.Delta;
    Real const beta  = [beta = desc.beta]() noexcept { // avoid beta == 1
        constexpr Real eps  = 1e-5;
        Real const     diff = beta - 1;
        return beta + (std::abs(diff) < eps ? std::copysign(eps, diff) : 0);
    }();
    rs    = RejectionSampler{ Delta, beta };
    vth1  = std::sqrt(desc.beta1) * c * std::abs(desc.Oc) / desc.op;
    T2OT1 = desc.T2_T1;
    xd    = desc.Vd / vth1;
    //
    xth2_squared = T2OT1 / (1 + (1 - Delta) * beta);
    vth1_cubed   = vth1 * vth1 * vth1;
}

auto LossconeVDF::f0(Vector const &v) const noexcept -> Real
{
    // note that vel = {v1, v2, v3}/vth1
    //
    // f0(x1, x2, x3) = exp(-(x1 - xd)^2)/√π *
    // ((1 - Δ*β)*exp(-(x2^2 + x3^2)/xth2^2) - (1 - Δ)*exp(-(x2^2 + x3^2)/(β*xth2^2)))
    // -------------------------------------------------------------------------------
    //                           (π * xth2^2 * (1 - β))
    //
    using std::exp;
    Real const f1 = [t = v.x - xd]() noexcept {
        return exp(-t * t) * M_2_SQRTPI * .5;
    }();
    Real const f2 = [D = rs.Delta, b = rs.beta, t2 = (v.y * v.y + v.z * v.z) / xth2_squared,
                     de = M_PI * xth2_squared * (1 - rs.beta)]() noexcept {
        return ((1 - D * b) * exp(-t2) - (1 - D) * exp(-t2 / b)) / de;
    }();
    return f1 * f2;
}

auto LossconeVDF::impl_emit() const -> Particle
{
    Particle ptl = load();

    // rescale
    //
    ptl.vel *= vth1;
    ptl.pos_x *= domain_extent.len;
    ptl.pos_x += domain_extent.loc;

    // delta-f parameters
    //
    ptl.psd.f = f0(ptl);
    // ptl.fOg = ptl.f/ptl.g0(ptl);
    static_assert(Particle::PSD::fOg == 1.0, "f and g should be identical");

    return ptl;
}
auto LossconeVDF::load() const -> Particle
{
    // position
    //
    Real const pos_x = bit_reversed<2>(); // [0, 1]

    // velocity in field-aligned frame (Hu et al., 2010, doi:10.1029/2009JA015158)
    //
    Real const phi1 = bit_reversed<3>() * 2 * M_PI;                               // [0, 2pi]
    Real const v1   = std::sqrt(-std::log(uniform_real<100>())) * std::sin(phi1); // v_para
    //
    auto const [v2, v3] = [v2 = rs.sample() * std::sqrt(xth2_squared)]() noexcept {
        Real const phi2 = bit_reversed<5>() * 2 * M_PI; // [0, 2pi]
        // {in-plane v_perp, out-of-plane v_perp}
        return std::make_pair(std::cos(phi2) * v2, std::sin(phi2) * v2);
    }();

    // velocity in Cartesian frame
    //
    Vector const vel = geomtr.fac2cart({ v1 + xd, v2, v3 });

    return Particle{ vel, pos_x };
}

// MARK: - RejectionSampler
//
LossconeVDF::RejectionSampler::RejectionSampler(Real const Delta, Real const beta /*must not be 1*/)
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
auto LossconeVDF::RejectionSampler::fOg(const Real x) const noexcept -> Real
{
    using std::exp;
    Real const x2 = x * x;
    Real const f  = ((1 - Delta * beta) * exp(-x2) - (1 - Delta) * exp(-x2 / beta)) / (1 - beta);
    Real const g  = exp(-x2 / alpha) / alpha;
    return f / g; // ratio of the target distribution to proposed distribution
}
auto LossconeVDF::RejectionSampler::sample() const noexcept -> Real
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
