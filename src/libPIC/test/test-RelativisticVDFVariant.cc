/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "catch2/catch.hpp"

#include <PIC/RelativisticVDFVariant.h>
#include <cmath>

using Particle = RelativisticParticle;

namespace {
bool operator==(Vector const &a, Vector const &b)
{
    return a.x == b.x && a.y == b.y && a.z == b.z;
}
bool operator==(CurviCoord const &a, CurviCoord const &b)
{
    return a.q1 == b.q1;
}
} // namespace
using ::operator==;

TEST_CASE("Test libPIC::RelativisticVDFVariant::RelativisticTestParticleVDF", "[libPIC::RelativisticVDFVariant::RelativisticTestParticleVDF]")
{
    Real const O0 = 1., op = 4 * O0, c = op;
    Real const xi = 0, D1 = 1;
    long const q1min = -7, q1max = 15;
    auto const geo   = Geometry{ xi, D1, O0 };
    auto const Nptls = 2;
    auto const desc  = TestParticleDesc<Nptls>(
        { -O0, op },
        { Vector{ 1, 2, 3 }, { 3, 4, 5 } },
        { CurviCoord{ q1min }, CurviCoord{ q1max } });
    auto const vdf = RelativisticVDFVariant::make(desc, geo, Range{ q1min, q1max - q1min }, c);

    CHECK(serialize(desc) == serialize(vdf.plasma_desc()));

    CHECK(vdf.Nrefcell_div_Ntotal() == Approx{ 1.0 }.epsilon(1e-10));

    for (long q1 = q1min; q1 <= q1max; ++q1) {
        CurviCoord const pos{ Real(q1) };

        auto const n0_ref = 0;
        auto const n0     = vdf.n0(pos);
        CHECK(*n0 == Approx{ n0_ref }.margin(1e-10));

        auto const nV0 = geo.cart_to_fac(vdf.nV0(pos), pos);
        CHECK(nV0.x == Approx{ 0 }.margin(1e-10));
        CHECK(nV0.y == Approx{ 0 }.margin(1e-10));
        CHECK(nV0.z == Approx{ 0 }.margin(1e-10));

        auto const nuv0 = geo.cart_to_fac(vdf.nuv0(pos), pos);
        REQUIRE(nuv0.tt == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ts.x == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ts.y == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ts.z == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.xx == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.yy == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.zz == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.xy == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.yz == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.zx == Approx{ 0 }.margin(1e-10));
    }

    // sampling
    auto const n_samples = Nptls + 1;
    auto const particles = vdf.emit(n_samples);

    for (unsigned i = 0; i < Nptls; ++i) {
        auto const &ptl = particles[i];

        REQUIRE(ptl.psd.weight == 1);
        REQUIRE(ptl.psd.real_f == -1);
        REQUIRE(ptl.psd.marker == -1);
        REQUIRE(ptl.gamma == Approx{ std::sqrt(1 + dot(ptl.g_vel, ptl.g_vel) / (c * c)) }.epsilon(1e-10));

        REQUIRE(ptl.g_vel == desc.vel[i]);
        REQUIRE(ptl.pos == desc.pos[i]);
    }
    {
        auto const &ptl = particles.back();

        REQUIRE(std::isnan(ptl.g_vel.x));
        REQUIRE(std::isnan(ptl.g_vel.y));
        REQUIRE(std::isnan(ptl.g_vel.z));
        REQUIRE(std::isnan(ptl.pos.q1));
        REQUIRE(std::isnan(ptl.gamma));
    }
}

TEST_CASE("Test libPIC::RelativisticVDFVariant::MaxwellianVDF", "[libPIC::RelativisticVDFVariant::MaxwellianVDF]")
{
    Real const O0 = 1., op = 4 * O0, c = op, beta1_eq = .1, T2OT1_eq = 5.35;
    Real const xi = 0, D1 = 1;
    long const q1min = -7, q1max = 15;
    auto const geo  = Geometry{ xi, D1, O0 };
    auto const desc = BiMaxPlasmaDesc({ { -O0, op }, 10, ShapeOrder::CIC }, beta1_eq, T2OT1_eq);
    auto const vdf  = RelativisticMaxwellianVDF(desc, geo, { q1min, q1max - q1min }, c);

    CHECK(serialize(desc) == serialize(vdf.plasma_desc()));

    // check equilibrium macro variables
    CHECK(vdf.Nrefcell_div_Ntotal() == Approx{ 1.0 / (q1max - q1min) }.epsilon(1e-10));

    for (long q1 = q1min; q1 <= q1max; ++q1) {
        CurviCoord const pos{ Real(q1) };
        auto const       eta = 1;

        auto const n0_ref = eta;
        auto const n0     = vdf.n0(pos);
        REQUIRE(*n0 == Approx{ n0_ref }.epsilon(1e-10));

        auto const nV0 = geo.cart_to_fac(vdf.nV0(pos), pos);
        REQUIRE(nV0.x == Approx{ 0 }.margin(1e-10));
        REQUIRE(nV0.y == Approx{ 0 }.margin(1e-10));
        REQUIRE(nV0.z == Approx{ 0 }.margin(1e-10));

        auto const nuv0 = geo.cart_to_fac(vdf.nuv0(pos), pos);
        REQUIRE(nuv0.tt == Approx{ *n0 * (c * c + .5 * beta1_eq * (.5 + T2OT1_eq * eta)) }.epsilon(1e-10));
        REQUIRE(nuv0.ts.x == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ts.y == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ts.z == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.xx == Approx{ *n0 * beta1_eq / 2 * (1 - beta1_eq / (c * c) * (3. / 4 + .5 * T2OT1_eq * eta)) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.yy == Approx{ *n0 * T2OT1_eq * eta * beta1_eq / 2 * (1 - beta1_eq / (c * c) * (1. / 4 + T2OT1_eq * eta)) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.zz == Approx{ *n0 * T2OT1_eq * eta * beta1_eq / 2 * (1 - beta1_eq / (c * c) * (1. / 4 + T2OT1_eq * eta)) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.xy == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.yz == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.zx == Approx{ 0 }.margin(1e-10));
    }

    // sampling
    auto const n_samples = 100U;
    auto const particles = vdf.emit(n_samples);

    std::for_each_n(begin(particles), n_samples, [c](Particle const &ptl) {
        REQUIRE(ptl.psd.weight == 1);
        REQUIRE(ptl.gamma == Approx{ std::sqrt(1 + dot(ptl.g_vel, ptl.g_vel) / (c * c)) }.epsilon(1e-10));
    });
}

TEST_CASE("Test libPIC::RelativisticVDFVariant::LossconeVDF", "[libPIC::RelativisticVDFVariant::LossconeVDF]")
{
    Real const O0 = 1., op = 4 * O0, c = op, beta1_eq = .1, T2OT1_eq = 5.35, beta_eq = .9;
    Real const xi = 0, D1 = 1;
    long const q1min = -7, q1max = 15;
    auto const geo     = Geometry{ xi, D1, O0 };
    auto const kinetic = KineticPlasmaDesc{ { -O0, op }, 10, ShapeOrder::CIC };
    auto const desc    = LossconePlasmaDesc(kinetic, beta1_eq, T2OT1_eq / (1 + beta_eq), beta_eq);
    auto const vdf     = RelativisticLossconeVDF(desc, geo, { q1min, q1max - q1min }, c);

    CHECK(serialize(desc) == serialize(vdf.plasma_desc()));

    // check equilibrium macro variables
    CHECK(vdf.Nrefcell_div_Ntotal() == Approx{ 1.0 / (q1max - q1min) }.epsilon(1e-10));

    for (long q1 = q1min; q1 <= q1max; ++q1) {
        CurviCoord const pos{ Real(q1) };
        auto const       eta               = 1;
        auto const       eta_b             = 1;
        auto const       beta              = beta_eq * eta_b / eta;
        auto const       vth_ratio_squared = T2OT1_eq / (1 + beta_eq) * eta;

        auto const n0_ref = (eta - beta_eq * eta_b) / (1 - beta_eq);
        auto const n0     = vdf.n0(pos);
        REQUIRE(*n0 == Approx{ n0_ref }.epsilon(1e-10));

        auto const nV0 = geo.cart_to_fac(vdf.nV0(pos), pos);
        REQUIRE(nV0.x == Approx{ 0 }.margin(1e-10));
        REQUIRE(nV0.y == Approx{ 0 }.margin(1e-10));
        REQUIRE(nV0.z == Approx{ 0 }.margin(1e-10));

        auto const nuv0 = geo.cart_to_fac(vdf.nuv0(pos), pos);
        REQUIRE(nuv0.tt == Approx{ *n0 * (c * c + .5 * beta1_eq * (.5 + vth_ratio_squared * (1 + beta))) }.epsilon(1e-10));
        REQUIRE(nuv0.ts.x == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ts.y == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ts.z == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.xx == Approx{ *n0 * beta1_eq / 2 * (1 - beta1_eq / (c * c) * (3. / 4 + .5 * vth_ratio_squared * (1 + beta))) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.yy == Approx{ *n0 * vth_ratio_squared * (1 + beta) * beta1_eq / 2 * (1 - beta1_eq / (c * c) * (1. / 4 + vth_ratio_squared * (1 + (1 + beta) * beta) / (1 + beta))) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.zz == Approx{ *n0 * vth_ratio_squared * (1 + beta) * beta1_eq / 2 * (1 - beta1_eq / (c * c) * (1. / 4 + vth_ratio_squared * (1 + (1 + beta) * beta) / (1 + beta))) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.xy == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.yz == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.zx == Approx{ 0 }.margin(1e-10));
    }

    // sampling
    auto const n_samples = 100U;
    auto const particles = vdf.emit(n_samples);

    std::for_each_n(begin(particles), n_samples, [c](Particle const &ptl) {
        REQUIRE(ptl.psd.weight == 1);
        REQUIRE(ptl.gamma == Approx{ std::sqrt(1 + dot(ptl.g_vel, ptl.g_vel) / (c * c)) }.epsilon(1e-10));
    });
}

TEST_CASE("Test libPIC::VDFVariant::RelativisticPartialShellVDF", "[libPIC::VDFVariant::RelativisticPartialShellVDF]")
{
    Real const O0 = 1, op = 4 * O0, c = op, beta = 0.1;
    Real const xi = 0, D1 = 1;
    long const q1min = -7, q1max = 15;
    auto const geo  = Geometry{ xi, D1, O0 };
    auto const desc = PartialShellPlasmaDesc({ { -O0, op }, 10, ShapeOrder::CIC }, beta);
    auto const vdf  = RelativisticVDFVariant::make(desc, geo, Range{ q1min, q1max - q1min }, c);

    CHECK(serialize(static_cast<KineticPlasmaDesc const &>(desc)) == serialize(vdf.plasma_desc()));

    // check equilibrium macro variables
    CHECK(vdf.Nrefcell_div_Ntotal() == Approx{ 1.0 / (q1max - q1min) }.epsilon(1e-10));

    for (long q1 = q1min; q1 <= q1max; ++q1) {
        CurviCoord const pos{ Real(q1) };
        auto const       eta = 1;

        auto const n0_ref = eta;
        auto const n0     = vdf.n0(pos);
        CHECK(*n0 == Approx{ n0_ref }.epsilon(1e-10));

        auto const nV0 = geo.cart_to_fac(vdf.nV0(pos), pos);
        CHECK(nV0.x == Approx{ 0 }.margin(1e-10));
        CHECK(nV0.y == Approx{ 0 }.margin(1e-10));
        CHECK(nV0.z == Approx{ 0 }.margin(1e-10));

        auto const nuv0 = geo.cart_to_fac(vdf.nuv0(pos), pos);
        auto const xs   = desc.vs / std::sqrt(beta);
        auto const xs2  = xs * xs;
        auto const Ab   = .5 * (xs * std::exp(-xs2) + 2 / M_2_SQRTPI * (.5 + xs2) * std::erfc(-xs));
        auto const T    = .5 * desc.beta / Ab * (xs * (2.5 + xs2) * std::exp(-xs2) + 2 / M_2_SQRTPI * (.75 + xs2 * (3 + xs2)) * std::erfc(-xs));
        auto const TT   = .5 * desc.beta * desc.beta / Ab * (xs * (33. / 4 + xs2 * (7 + xs2)) * std::exp(-xs2) + 2 / M_2_SQRTPI * (15. / 8 + xs2 * (45. / 4 + 15. / 2 * xs2 + xs2 * xs2)) * std::erfc(-xs));
        REQUIRE(nuv0.tt == Approx{ *n0 * (c * c + .5 * T) }.epsilon(1e-10));
        REQUIRE(nuv0.ts.x == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ts.y == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ts.z == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.xx == Approx{ *n0 / Real(3 + desc.zeta) * (T - .5 * TT / (c * c)) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.yy == Approx{ *n0 * (1 + .5 * desc.zeta) / Real(3 + desc.zeta) * (T - .5 * TT / (c * c)) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.zz == Approx{ *n0 * (1 + .5 * desc.zeta) / Real(3 + desc.zeta) * (T - .5 * TT / (c * c)) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.xy == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.yz == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.zx == Approx{ 0 }.margin(1e-10));
    }

    // sampling
    auto const n_samples = 100U;
    auto const particles = vdf.emit(n_samples);

    std::for_each_n(begin(particles), n_samples, [c](Particle const &ptl) {
        REQUIRE(ptl.psd.weight == 1);
        REQUIRE(ptl.gamma == Approx{ std::sqrt(1 + dot(ptl.g_vel, ptl.g_vel) / (c * c)) }.epsilon(1e-10));
    });
}
