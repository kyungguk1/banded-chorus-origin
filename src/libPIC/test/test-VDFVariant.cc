/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "catch2/catch.hpp"

#include <PIC/VDFVariant.h>
#include <cmath>

TEST_CASE("Test libPIC::VDFVariant::MaxwellianVDF", "[libPIC::VDFVariant::MaxwellianVDF]")
{
    Real const O0 = 1., op = 4 * O0, c = op, beta1_eq = .1, T2OT1_eq = 5.35;
    Real const xi = 0, D1 = 1;
    long const q1min = -7, q1max = 15;
    auto const geo  = Geometry{ xi, D1, O0 };
    auto const desc = BiMaxPlasmaDesc({ { -O0, op }, 10, ShapeOrder::CIC }, beta1_eq, T2OT1_eq);
    auto const vdf  = VDFVariant::make(desc, geo, Range{ q1min, q1max - q1min }, c);

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

        auto const nvv0 = geo.cart_to_fac(vdf.nvv0(pos), pos);
        CHECK(nvv0.xx == Approx{ desc.beta1 / 2 * eta }.epsilon(1e-10));
        CHECK(nvv0.yy == Approx{ desc.beta1 / 2 * desc.T2_T1 * eta * eta }.epsilon(1e-10));
        CHECK(nvv0.zz == Approx{ desc.beta1 / 2 * desc.T2_T1 * eta * eta }.epsilon(1e-10));
        CHECK(nvv0.xy == Approx{ 0 }.margin(1e-10));
        CHECK(nvv0.yz == Approx{ 0 }.margin(1e-10));
        CHECK(nvv0.zx == Approx{ 0 }.margin(1e-10));
    }

    // sampling
    auto const n_samples = 100U;
    auto const particles = vdf.emit(n_samples);

    std::for_each_n(begin(particles), n_samples, [](Particle const &ptl) {
        REQUIRE(ptl.psd.weight == 1);
        REQUIRE(ptl.psd.real_f == -1);
        REQUIRE(ptl.psd.marker == -1);
    });
}

TEST_CASE("Test libPIC::VDFVariant::LossconeVDF::Loss", "[libPIC::VDFVariant::LossconeVDF::Loss]")
{
    Real const O0 = 1., op = 4 * O0, c = op, beta1_eq = .1, T2OT1_eq = 5.35, beta_eq = .9;
    Real const xi = 0, D1 = 1;
    long const q1min = -7, q1max = 15;
    auto const geo     = Geometry{ xi, D1, O0 };
    auto const kinetic = KineticPlasmaDesc{ { -O0, op }, 10, ShapeOrder::CIC };
    auto const desc    = LossconePlasmaDesc(kinetic, beta1_eq, T2OT1_eq / (1 + beta_eq), beta_eq);
    auto const vdf     = VDFVariant::make(desc, geo, Range{ q1min, q1max - q1min }, c);

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

        auto const nvv0 = geo.cart_to_fac(vdf.nvv0(pos), pos);
        CHECK(nvv0.xx == Approx{ desc.beta1 / 2 * eta }.epsilon(1e-10));
        CHECK(nvv0.yy == Approx{ desc.beta1 / 2 * desc.T2_T1 * eta * eta }.epsilon(1e-10));
        CHECK(nvv0.zz == Approx{ desc.beta1 / 2 * desc.T2_T1 * eta * eta }.epsilon(1e-10));
        CHECK(nvv0.xy == Approx{ 0 }.margin(1e-10));
        CHECK(nvv0.yz == Approx{ 0 }.margin(1e-10));
        CHECK(nvv0.zx == Approx{ 0 }.margin(1e-10));
    }

    // sampling
    auto const n_samples = 100U;
    auto const particles = vdf.emit(n_samples);

    std::for_each_n(begin(particles), n_samples, [](Particle const &ptl) {
        REQUIRE(ptl.psd.weight == 1);
        REQUIRE(ptl.psd.real_f == -1);
        REQUIRE(ptl.psd.marker == -1);
    });
}
