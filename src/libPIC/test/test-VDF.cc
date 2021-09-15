/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "catch2/catch.hpp"

#include <PIC/LossconeVDF.h>
#include <PIC/MaxwellianVDF.h>
#include <algorithm>
#include <cmath>

TEST_CASE("Test libPIC::MaxwellianVDF::full_f", "[libPIC::MaxwellianVDF::full_f]")
{
    auto const O0 = 1., theta = 30. * M_PI / 180, op = 4 * O0, c = op;
    auto const extent = Range{ 0, 10 } - 1;
    auto const geo    = Geometry{ O0, theta };
    auto const desc   = BiMaxPlasmaDesc({ { -O0, op }, 10, ShapeOrder::CIC }, .1, 2, -1);
    auto const vdf    = MaxwellianVDF(desc, geo, extent, c);

    CHECK(serialize(desc) == serialize(vdf.plasma_desc()));

    // check equilibrium macro variables
    auto const n0 = vdf.n0(0);
    CHECK(*n0 == Approx{ 1 }.epsilon(1e-10));

    auto const nV0 = geo.cart2fac(vdf.nV0(0));
    CHECK(nV0.x == Approx{ desc.Vd }.epsilon(1e-10));
    CHECK(nV0.y == Approx{ 0 }.margin(1e-10));
    CHECK(nV0.z == Approx{ 0 }.margin(1e-10));

    auto const nvv0 = geo.cart2fac(vdf.nvv0(0));
    CHECK(nvv0.xx == Approx{ desc.beta1 / 2 + Real{ n0 } * desc.Vd * desc.Vd }.epsilon(1e-10));
    CHECK(nvv0.yy == Approx{ desc.beta1 / 2 * desc.T2_T1 }.epsilon(1e-10));
    CHECK(nvv0.zz == Approx{ desc.beta1 / 2 * desc.T2_T1 }.epsilon(1e-10));
    CHECK(nvv0.xy == Approx{ 0 }.margin(1e-10));
    CHECK(nvv0.yz == Approx{ 0 }.margin(1e-10));
    CHECK(nvv0.zx == Approx{ 0 }.margin(1e-10));

    // sampling
    auto const n_samples = 100000U;
    auto const particles = vdf.emit(n_samples);

    auto pos_mom1 = Real{};
    auto pos_mom2 = Real{};
    auto vel_mom1 = Vector{};
    auto vel_mom2 = Tensor{};
    for (auto const &ptl : particles) {
        pos_mom1 += ptl.pos_x / n_samples;
        pos_mom2 += ptl.pos_x * ptl.pos_x / n_samples;
        vel_mom1 += ptl.vel / n_samples;
        Tensor vv;
        vv.lo() = vv.hi() = ptl.vel;
        vv.lo() *= ptl.vel;
        vv.hi() *= { ptl.vel.y, ptl.vel.z, ptl.vel.x };
        vel_mom2 += vv / n_samples;
    }
    vel_mom1 = geo.cart2fac(vel_mom1);
    vel_mom2 = geo.cart2fac(vel_mom2);

    auto const X_exact = extent.mean();
    CHECK(std::abs(pos_mom1 - X_exact) < X_exact * 1e-4);
    auto const xx_exact
        = (extent.min() * extent.min() + extent.min() * extent.max() + extent.max() * extent.max()) / 3.;
    CHECK(std::abs(pos_mom2 - xx_exact) < xx_exact * 1e-4);

    CHECK(vel_mom1.x == Approx{ nV0.x }.epsilon(1e-3));
    CHECK(vel_mom1.y == Approx{ 0 }.margin(1e-3));
    CHECK(vel_mom1.z == Approx{ 0 }.margin(1e-3));

    CHECK(vel_mom2.xx == Approx{ nvv0.xx }.epsilon(1e-3));
    CHECK(vel_mom2.yy == Approx{ nvv0.yy }.epsilon(4e-3));
    CHECK(vel_mom2.zz == Approx{ nvv0.zz }.epsilon(5e-3));
    CHECK(vel_mom2.xy == Approx{ nvv0.xy }.margin(1e-3));
    CHECK(vel_mom2.yz == Approx{ nvv0.yz }.margin(1e-3));
    CHECK(vel_mom2.zx == Approx{ nvv0.zx }.margin(1e-3));

    static_assert(n_samples > 10);
    std::for_each_n(begin(particles), 10, [](auto const &ptl) {
        REQUIRE(ptl.psd.weight == 1);
        REQUIRE(ptl.psd.real_f == -1);
        REQUIRE(ptl.psd.marker == -1);
    });
}

TEST_CASE("Test libPIC::MaxwellianVDF::delta_f", "[libPIC::MaxwellianVDF::delta_f]")
{
    auto const O0 = 1., theta = 30. * M_PI / 180, op = 4 * O0, c = op;
    auto const extent  = Range{ 0, 10 } - 1;
    auto const geo     = Geometry{ O0, theta };
    auto const kinetic = KineticPlasmaDesc{ { -O0, op }, 10, ShapeOrder::CIC, ParticleScheme::delta_f, .001, 2.1 };
    auto const desc    = BiMaxPlasmaDesc(kinetic, .1, 2, -1);
    auto const vdf     = MaxwellianVDF(desc, geo, extent, c);

    auto const g_desc = BiMaxPlasmaDesc({ { -O0, op }, 10, ShapeOrder::CIC }, .1 * desc.marker_temp_ratio, 2, -1);
    auto const g_vdf  = MaxwellianVDF(g_desc, geo, extent, c);

    CHECK(serialize(desc) == serialize(vdf.plasma_desc()));

    // check equilibrium macro variables
    auto const n0 = vdf.n0(0);
    CHECK(*n0 == Approx{ 1 }.epsilon(1e-10));

    auto const nV0 = geo.cart2fac(vdf.nV0(0));
    CHECK(nV0.x == Approx{ desc.Vd }.epsilon(1e-10));
    CHECK(nV0.y == Approx{ 0 }.margin(1e-10));
    CHECK(nV0.z == Approx{ 0 }.margin(1e-10));

    auto const nvv0 = geo.cart2fac(vdf.nvv0(0));
    CHECK(nvv0.xx == Approx{ desc.beta1 / 2 + Real{ n0 } * desc.Vd * desc.Vd }.epsilon(1e-10));
    CHECK(nvv0.yy == Approx{ desc.beta1 / 2 * desc.T2_T1 }.epsilon(1e-10));
    CHECK(nvv0.zz == Approx{ desc.beta1 / 2 * desc.T2_T1 }.epsilon(1e-10));
    CHECK(nvv0.xy == Approx{ 0 }.margin(1e-10));
    CHECK(nvv0.yz == Approx{ 0 }.margin(1e-10));
    CHECK(nvv0.zx == Approx{ 0 }.margin(1e-10));

    // sampling
    auto const n_samples = 100000U;
    auto const particles = vdf.emit(n_samples);

    auto pos_mom1 = Real{};
    auto pos_mom2 = Real{};
    auto vel_mom1 = Vector{};
    auto vel_mom2 = Tensor{};
    for (auto const &ptl : particles) {
        pos_mom1 += ptl.pos_x / n_samples;
        pos_mom2 += ptl.pos_x * ptl.pos_x / n_samples;
        vel_mom1 += ptl.vel / n_samples;
        Tensor vv;
        vv.lo() = vv.hi() = ptl.vel;
        vv.lo() *= ptl.vel;
        vv.hi() *= { ptl.vel.y, ptl.vel.z, ptl.vel.x };
        vel_mom2 += vv / n_samples;
    }
    vel_mom1 = geo.cart2fac(vel_mom1);
    vel_mom2 = geo.cart2fac(vel_mom2);

    auto const X_exact = extent.mean();
    CHECK(std::abs(pos_mom1 - X_exact) < X_exact * 1e-4);
    auto const xx_exact
        = (extent.min() * extent.min() + extent.min() * extent.max() + extent.max() * extent.max()) / 3.;
    CHECK(std::abs(pos_mom2 - xx_exact) < xx_exact * 1e-4);

    auto const g_nV0 = geo.cart2fac(g_vdf.nV0(0));
    CHECK(vel_mom1.x == Approx{ g_nV0.x }.epsilon(1e-3));
    CHECK(vel_mom1.y == Approx{ 0 }.margin(1e-3));
    CHECK(vel_mom1.z == Approx{ 0 }.margin(1e-3));

    auto const g_nvv0 = geo.cart2fac(g_vdf.nvv0(0));
    CHECK(vel_mom2.xx == Approx{ g_nvv0.xx }.epsilon(2e-3));
    CHECK(vel_mom2.yy == Approx{ g_nvv0.yy }.epsilon(1e-3));
    CHECK(vel_mom2.zz == Approx{ g_nvv0.zz }.epsilon(2e-3));
    CHECK(vel_mom2.xy == Approx{ g_nvv0.xy }.margin(1e-3));
    CHECK(vel_mom2.yz == Approx{ g_nvv0.yz }.margin(1e-3));
    CHECK(vel_mom2.zx == Approx{ g_nvv0.zx }.margin(1e-3));

    static_assert(n_samples > 10);
    std::for_each_n(begin(particles), 10, [&](auto const &ptl) {
        REQUIRE(ptl.psd.weight == desc.initial_weight);
        REQUIRE(ptl.psd.marker == g_vdf.f0(ptl));
        REQUIRE(ptl.psd.real_f == Approx{ vdf.f0(ptl) + ptl.psd.weight * ptl.psd.marker }.epsilon(1e-10));
        REQUIRE(ptl.psd.weight == Approx{ vdf.weight(ptl) }.epsilon(1e-10));
    });
}

TEST_CASE("Test libPIC::LossconeVDF::BiMax::full_f", "[libPIC::LossconeVDF::BiMax::full_f]")
{
    auto const O0 = 1., theta = 30. * M_PI / 180, op = 4 * O0, c = op;
    auto const extent = Range{ 0, 10 } - 1;
    auto const geo    = Geometry{ O0, theta };
    auto const desc   = BiMaxPlasmaDesc({ { -O0, op }, 10, ShapeOrder::CIC }, .1, 2, -1);
    auto const vdf    = LossconeVDF(LossconePlasmaDesc{ desc }, geo, extent, c);

    CHECK(serialize(LossconePlasmaDesc{ desc }) == serialize(vdf.plasma_desc()));

    // check equilibrium macro variables
    auto const n0 = vdf.n0(0);
    CHECK(*n0 == Approx{ 1 }.epsilon(1e-10));

    auto const nV0 = geo.cart2fac(vdf.nV0(0));
    CHECK(nV0.x == Approx{ desc.Vd }.epsilon(1e-10));
    CHECK(nV0.y == Approx{ 0 }.margin(1e-10));
    CHECK(nV0.z == Approx{ 0 }.margin(1e-10));

    auto const nvv0 = geo.cart2fac(vdf.nvv0(0));
    CHECK(nvv0.xx == Approx{ desc.beta1 / 2 + Real{ n0 } * desc.Vd * desc.Vd }.epsilon(1e-10));
    CHECK(nvv0.yy == Approx{ desc.beta1 / 2 * desc.T2_T1 }.epsilon(1e-10));
    CHECK(nvv0.zz == Approx{ desc.beta1 / 2 * desc.T2_T1 }.epsilon(1e-10));
    CHECK(nvv0.xy == Approx{ 0 }.margin(1e-10));
    CHECK(nvv0.yz == Approx{ 0 }.margin(1e-10));
    CHECK(nvv0.zx == Approx{ 0 }.margin(1e-10));

    // sampling
    auto const n_samples = 100000U;
    auto const particles = vdf.emit(n_samples);

    auto pos_mom1 = Real{};
    auto pos_mom2 = Real{};
    auto vel_mom1 = Vector{};
    auto vel_mom2 = Tensor{};
    for (auto const &ptl : particles) {
        pos_mom1 += ptl.pos_x / n_samples;
        pos_mom2 += ptl.pos_x * ptl.pos_x / n_samples;
        vel_mom1 += ptl.vel / n_samples;
        Tensor vv;
        vv.lo() = vv.hi() = ptl.vel;
        vv.lo() *= ptl.vel;
        vv.hi() *= { ptl.vel.y, ptl.vel.z, ptl.vel.x };
        vel_mom2 += vv / n_samples;
    }
    vel_mom1 = geo.cart2fac(vel_mom1);
    vel_mom2 = geo.cart2fac(vel_mom2);

    auto const X_exact = extent.mean();
    CHECK(std::abs(pos_mom1 - X_exact) < X_exact * 1e-4);
    auto const xx_exact
        = (extent.min() * extent.min() + extent.min() * extent.max() + extent.max() * extent.max()) / 3.;
    CHECK(std::abs(pos_mom2 - xx_exact) < xx_exact * 1e-4);

    CHECK(vel_mom1.x == Approx{ nV0.x }.epsilon(1e-3));
    CHECK(vel_mom1.y == Approx{ 0 }.margin(1e-3));
    CHECK(vel_mom1.z == Approx{ 0 }.margin(1e-3));

    CHECK(vel_mom2.xx == Approx{ nvv0.xx }.epsilon(1e-3));
    CHECK(vel_mom2.yy == Approx{ nvv0.yy }.epsilon(4e-3));
    CHECK(vel_mom2.zz == Approx{ nvv0.zz }.epsilon(5e-3));
    CHECK(vel_mom2.xy == Approx{ nvv0.xy }.margin(1e-3));
    CHECK(vel_mom2.yz == Approx{ nvv0.yz }.margin(1e-3));
    CHECK(vel_mom2.zx == Approx{ nvv0.zx }.margin(1e-3));

    static_assert(n_samples > 10);
    std::for_each_n(begin(particles), 10, [](auto const &ptl) {
        REQUIRE(ptl.psd.weight == 1);
        REQUIRE(ptl.psd.real_f == -1);
        REQUIRE(ptl.psd.marker == -1);
    });
}

TEST_CASE("Test libPIC::LossconeVDF::BiMax::delta_f", "[libPIC::LossconeVDF::BiMax::delta_f]")
{
    auto const O0 = 1., theta = 30. * M_PI / 180, op = 4 * O0, c = op;
    auto const extent  = Range{ 0, 10 } - 1;
    auto const geo     = Geometry{ O0, theta };
    auto const kinetic = KineticPlasmaDesc{ { -O0, op }, 10, ShapeOrder::CIC, ParticleScheme::delta_f, .001, 2.1 };
    auto const desc    = BiMaxPlasmaDesc(kinetic, .1, 2, -1);
    auto const vdf     = LossconeVDF(LossconePlasmaDesc{ desc }, geo, extent, c);

    auto const g_desc = BiMaxPlasmaDesc({ { -O0, op }, 10, ShapeOrder::CIC }, .1 * desc.marker_temp_ratio, 2, -1);
    auto const g_vdf  = LossconeVDF(LossconePlasmaDesc{ g_desc }, geo, extent, c);

    CHECK(serialize(LossconePlasmaDesc{ desc }) == serialize(vdf.plasma_desc()));

    // check equilibrium macro variables
    auto const n0 = vdf.n0(0);
    CHECK(*n0 == Approx{ 1 }.epsilon(1e-10));

    auto const nV0 = geo.cart2fac(vdf.nV0(0));
    CHECK(nV0.x == Approx{ desc.Vd }.epsilon(1e-10));
    CHECK(nV0.y == Approx{ 0 }.margin(1e-10));
    CHECK(nV0.z == Approx{ 0 }.margin(1e-10));

    auto const nvv0 = geo.cart2fac(vdf.nvv0(0));
    CHECK(nvv0.xx == Approx{ desc.beta1 / 2 + Real{ n0 } * desc.Vd * desc.Vd }.epsilon(1e-10));
    CHECK(nvv0.yy == Approx{ desc.beta1 / 2 * desc.T2_T1 }.epsilon(1e-10));
    CHECK(nvv0.zz == Approx{ desc.beta1 / 2 * desc.T2_T1 }.epsilon(1e-10));
    CHECK(nvv0.xy == Approx{ 0 }.margin(1e-10));
    CHECK(nvv0.yz == Approx{ 0 }.margin(1e-10));
    CHECK(nvv0.zx == Approx{ 0 }.margin(1e-10));

    // sampling
    auto const n_samples = 100000U;
    auto const particles = vdf.emit(n_samples);

    auto pos_mom1 = Real{};
    auto pos_mom2 = Real{};
    auto vel_mom1 = Vector{};
    auto vel_mom2 = Tensor{};
    for (auto const &ptl : particles) {
        pos_mom1 += ptl.pos_x / n_samples;
        pos_mom2 += ptl.pos_x * ptl.pos_x / n_samples;
        vel_mom1 += ptl.vel / n_samples;
        Tensor vv;
        vv.lo() = vv.hi() = ptl.vel;
        vv.lo() *= ptl.vel;
        vv.hi() *= { ptl.vel.y, ptl.vel.z, ptl.vel.x };
        vel_mom2 += vv / n_samples;
    }
    vel_mom1 = geo.cart2fac(vel_mom1);
    vel_mom2 = geo.cart2fac(vel_mom2);

    auto const X_exact = extent.mean();
    CHECK(std::abs(pos_mom1 - X_exact) < X_exact * 1e-4);
    auto const xx_exact
        = (extent.min() * extent.min() + extent.min() * extent.max() + extent.max() * extent.max()) / 3.;
    CHECK(std::abs(pos_mom2 - xx_exact) < xx_exact * 1e-4);

    auto const g_nV0 = geo.cart2fac(g_vdf.nV0(0));
    CHECK(vel_mom1.x == Approx{ g_nV0.x }.epsilon(1e-3));
    CHECK(vel_mom1.y == Approx{ 0 }.margin(1e-3));
    CHECK(vel_mom1.z == Approx{ 0 }.margin(1e-3));

    auto const g_nvv0 = geo.cart2fac(g_vdf.nvv0(0));
    CHECK(vel_mom2.xx == Approx{ g_nvv0.xx }.epsilon(2e-3));
    CHECK(vel_mom2.yy == Approx{ g_nvv0.yy }.epsilon(1e-3));
    CHECK(vel_mom2.zz == Approx{ g_nvv0.zz }.epsilon(2e-3));
    CHECK(vel_mom2.xy == Approx{ g_nvv0.xy }.margin(1e-3));
    CHECK(vel_mom2.yz == Approx{ g_nvv0.yz }.margin(1e-3));
    CHECK(vel_mom2.zx == Approx{ g_nvv0.zx }.margin(1e-3));

    static_assert(n_samples > 10);
    std::for_each_n(begin(particles), 10, [&](auto const &ptl) {
        REQUIRE(ptl.psd.weight == desc.initial_weight);
        REQUIRE(ptl.psd.marker == g_vdf.f0(ptl));
        REQUIRE(ptl.psd.real_f == Approx{ vdf.f0(ptl) + ptl.psd.weight * ptl.psd.marker }.epsilon(1e-10));
        REQUIRE(ptl.psd.weight == Approx{ vdf.weight(ptl) }.epsilon(1e-10));
    });
}

TEST_CASE("Test libPIC::LossconeVDF::Loss::full_f", "[libPIC::LossconeVDF::Loss::full_f]")
{
    auto const O0 = 1., theta = 30. * M_PI / 180, op = 4 * O0, c = op;
    auto const extent  = Range{ 0, 10 } - 1;
    auto const geo     = Geometry{ O0, theta };
    auto const kinetic = KineticPlasmaDesc{ { -O0, op }, 10, ShapeOrder::CIC };
    auto const desc    = LossconePlasmaDesc(kinetic, .1, 2, { .5, .9 }, 1.1);
    auto const vdf     = LossconeVDF(desc, geo, extent, c);

    CHECK(serialize(desc) == serialize(vdf.plasma_desc()));

    // check equilibrium macro variables
    auto const n0 = vdf.n0(0);
    CHECK(*n0 == Approx{ 1 }.epsilon(1e-10));

    auto const nV0 = geo.cart2fac(vdf.nV0(0));
    CHECK(nV0.x == Approx{ desc.Vd }.epsilon(1e-10));
    CHECK(nV0.y == Approx{ 0 }.margin(1e-10));
    CHECK(nV0.z == Approx{ 0 }.margin(1e-10));

    auto const nvv0 = geo.cart2fac(vdf.nvv0(0));
    CHECK(nvv0.xx == Approx{ desc.beta1 / 2 + Real{ n0 } * desc.Vd * desc.Vd }.epsilon(1e-10));
    CHECK(nvv0.yy == Approx{ desc.beta1 / 2 * desc.T2_T1 }.epsilon(1e-10));
    CHECK(nvv0.zz == Approx{ desc.beta1 / 2 * desc.T2_T1 }.epsilon(1e-10));
    CHECK(nvv0.xy == Approx{ 0 }.margin(1e-10));
    CHECK(nvv0.yz == Approx{ 0 }.margin(1e-10));
    CHECK(nvv0.zx == Approx{ 0 }.margin(1e-10));

    // sampling
    auto const n_samples = 100000U;
    auto const particles = vdf.emit(n_samples);

    auto pos_mom1 = Real{};
    auto pos_mom2 = Real{};
    auto vel_mom1 = Vector{};
    auto vel_mom2 = Tensor{};
    for (auto const &ptl : particles) {
        pos_mom1 += ptl.pos_x / n_samples;
        pos_mom2 += ptl.pos_x * ptl.pos_x / n_samples;
        vel_mom1 += ptl.vel / n_samples;
        Tensor vv;
        vv.lo() = vv.hi() = ptl.vel;
        vv.lo() *= ptl.vel;
        vv.hi() *= { ptl.vel.y, ptl.vel.z, ptl.vel.x };
        vel_mom2 += vv / n_samples;
    }
    vel_mom1 = geo.cart2fac(vel_mom1);
    vel_mom2 = geo.cart2fac(vel_mom2);

    auto const X_exact = extent.mean();
    CHECK(std::abs(pos_mom1 - X_exact) < X_exact * 1e-4);
    auto const xx_exact
        = (extent.min() * extent.min() + extent.min() * extent.max() + extent.max() * extent.max()) / 3.;
    CHECK(std::abs(pos_mom2 - xx_exact) < xx_exact * 1e-4);

    CHECK(vel_mom1.x == Approx{ nV0.x }.epsilon(1e-3));
    CHECK(vel_mom1.y == Approx{ 0 }.margin(1e-3));
    CHECK(vel_mom1.z == Approx{ 0 }.margin(1e-3));

    CHECK(vel_mom2.xx == Approx{ nvv0.xx }.epsilon(1e-3));
    CHECK(vel_mom2.yy == Approx{ nvv0.yy }.epsilon(3e-3));
    CHECK(vel_mom2.zz == Approx{ nvv0.zz }.epsilon(3e-3));
    CHECK(vel_mom2.xy == Approx{ nvv0.xy }.margin(1e-3));
    CHECK(vel_mom2.yz == Approx{ nvv0.yz }.margin(1e-3));
    CHECK(vel_mom2.zx == Approx{ nvv0.zx }.margin(1e-3));

    static_assert(n_samples > 10);
    std::for_each_n(begin(particles), 10, [](auto const &ptl) {
        REQUIRE(ptl.psd.weight == 1);
        REQUIRE(ptl.psd.real_f == -1);
        REQUIRE(ptl.psd.marker == -1);
    });
}

TEST_CASE("Test libPIC::LossconeVDF::Loss::delta_f", "[libPIC::LossconeVDF::Loss::delta_f]")
{
    auto const O0 = 1., theta = 30. * M_PI / 180, op = 4 * O0, c = op;
    auto const extent  = Range{ 0, 10 } - 1;
    auto const geo     = Geometry{ O0, theta };
    auto const kinetic = KineticPlasmaDesc{ { -O0, op }, 10, ShapeOrder::CIC, ParticleScheme::delta_f, .001, 2.1 };
    auto const desc    = LossconePlasmaDesc(kinetic, .1, 2, { .5, .9 }, 1.1);
    auto const vdf     = LossconeVDF(desc, geo, extent, c);

    auto const g_desc = LossconePlasmaDesc({ { -O0, op }, 10, ShapeOrder::CIC }, .1 * desc.marker_temp_ratio, 2, { .5, .9 }, 1.1);
    auto const g_vdf  = LossconeVDF(g_desc, geo, extent, c);

    CHECK(serialize(desc) == serialize(vdf.plasma_desc()));

    // check equilibrium macro variables
    auto const n0 = vdf.n0(0);
    CHECK(*n0 == Approx{ 1 }.epsilon(1e-10));

    auto const nV0 = geo.cart2fac(vdf.nV0(0));
    CHECK(nV0.x == Approx{ desc.Vd }.epsilon(1e-10));
    CHECK(nV0.y == Approx{ 0 }.margin(1e-10));
    CHECK(nV0.z == Approx{ 0 }.margin(1e-10));

    auto const nvv0 = geo.cart2fac(vdf.nvv0(0));
    CHECK(nvv0.xx == Approx{ desc.beta1 / 2 + Real{ n0 } * desc.Vd * desc.Vd }.epsilon(1e-10));
    CHECK(nvv0.yy == Approx{ desc.beta1 / 2 * desc.T2_T1 }.epsilon(1e-10));
    CHECK(nvv0.zz == Approx{ desc.beta1 / 2 * desc.T2_T1 }.epsilon(1e-10));
    CHECK(nvv0.xy == Approx{ 0 }.margin(1e-10));
    CHECK(nvv0.yz == Approx{ 0 }.margin(1e-10));
    CHECK(nvv0.zx == Approx{ 0 }.margin(1e-10));

    // sampling
    auto const n_samples = 100000U;
    auto const particles = vdf.emit(n_samples);

    auto pos_mom1 = Real{};
    auto pos_mom2 = Real{};
    auto vel_mom1 = Vector{};
    auto vel_mom2 = Tensor{};
    for (auto const &ptl : particles) {
        pos_mom1 += ptl.pos_x / n_samples;
        pos_mom2 += ptl.pos_x * ptl.pos_x / n_samples;
        vel_mom1 += ptl.vel / n_samples;
        Tensor vv;
        vv.lo() = vv.hi() = ptl.vel;
        vv.lo() *= ptl.vel;
        vv.hi() *= { ptl.vel.y, ptl.vel.z, ptl.vel.x };
        vel_mom2 += vv / n_samples;
    }
    vel_mom1 = geo.cart2fac(vel_mom1);
    vel_mom2 = geo.cart2fac(vel_mom2);

    auto const X_exact = extent.mean();
    CHECK(std::abs(pos_mom1 - X_exact) < X_exact * 1e-4);
    auto const xx_exact
        = (extent.min() * extent.min() + extent.min() * extent.max() + extent.max() * extent.max()) / 3.;
    CHECK(std::abs(pos_mom2 - xx_exact) < xx_exact * 1e-4);

    auto const g_nV0 = geo.cart2fac(g_vdf.nV0(0));
    CHECK(vel_mom1.x == Approx{ g_nV0.x }.epsilon(1e-3));
    CHECK(vel_mom1.y == Approx{ 0 }.margin(1e-3));
    CHECK(vel_mom1.z == Approx{ 0 }.margin(2e-3));

    auto const g_nvv0 = geo.cart2fac(g_vdf.nvv0(0));
    CHECK(vel_mom2.xx == Approx{ g_nvv0.xx }.epsilon(2e-3));
    CHECK(vel_mom2.yy == Approx{ g_nvv0.yy }.epsilon(3e-3));
    CHECK(vel_mom2.zz == Approx{ g_nvv0.zz }.epsilon(8e-3));
    CHECK(vel_mom2.xy == Approx{ g_nvv0.xy }.margin(1e-3));
    CHECK(vel_mom2.yz == Approx{ g_nvv0.yz }.margin(1e-3));
    CHECK(vel_mom2.zx == Approx{ g_nvv0.zx }.margin(1e-3));

    static_assert(n_samples > 10);
    std::for_each_n(begin(particles), 10, [&](auto const &ptl) {
        REQUIRE(ptl.psd.weight == desc.initial_weight);
        REQUIRE(ptl.psd.marker == g_vdf.f0(ptl));
        REQUIRE(ptl.psd.real_f == Approx{ vdf.f0(ptl) + ptl.psd.weight * ptl.psd.marker }.epsilon(1e-10));
        REQUIRE(ptl.psd.weight == Approx{ vdf.weight(ptl) }.epsilon(1e-10));
    });
}
