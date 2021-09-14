/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "catch2/catch.hpp"

#include <PIC/RelativisticVDFVariant.h>
#include <cmath>

TEST_CASE("Test libPIC::RelativisticVDFVariant::MaxwellianVDF", "[libPIC::RelativisticVDFVariant::MaxwellianVDF]")
{
    Real constexpr O0 = 1., c = 4., theta = 30. * M_PI / 180, op = c;
    auto constexpr extent = Range{ 0, 10 } - 1;
    auto constexpr desc   = BiMaxPlasmaDesc({ { -O0, op }, 10, ShapeOrder::CIC, ParticleScheme::delta_f, .01 },
                                            .1, 2, -1);
    auto const geo        = Geometry{ O0, theta };
    auto const vdf        = RelativisticVDFVariant::make(desc, geo, extent, c);

    CHECK(serialize(static_cast<KineticPlasmaDesc const &>(desc)) == serialize(vdf.plasma_desc()));

    auto const n0 = vdf.n0(0);
    CHECK(*n0 == Approx{ 1 }.epsilon(1e-10));

    auto const nV0 = geo.cart2fac(vdf.nV0(0));
    CHECK(nV0.x == Approx{ -0.999999998515408 }.epsilon(1e-8));
    CHECK(nV0.y == Approx{ 0 }.margin(1e-10));
    CHECK(nV0.z == Approx{ 0 }.margin(1e-10));

    auto const nuv0 = geo.cart2fac(vdf.nuv0(0));
    CHECK(nuv0.tt == Approx{ 16.65617028870888 }.epsilon(1e-3));
    CHECK(nuv0.ts.x == Approx{ -4.176016301152911 }.epsilon(1e-4));
    CHECK(nuv0.ts.y == Approx{ 0 }.margin(1e-10));
    CHECK(nuv0.ts.z == Approx{ 0 }.margin(1e-10));
    CHECK(nuv0.ss.xx == Approx{ 1.091899002238576 }.epsilon(1e-4));
    CHECK(nuv0.ss.yy == Approx{ 0.0955020983746306 }.epsilon(1e-3));
    CHECK(nuv0.ss.zz == Approx{ 0.0955020983746306 }.epsilon(1e-3));
    CHECK(nuv0.ss.xy == Approx{ 0 }.margin(1e-10));
    CHECK(nuv0.ss.yz == Approx{ 0 }.margin(1e-10));
    CHECK(nuv0.ss.zx == Approx{ 0 }.margin(1e-10));

    // sampling
    auto const n_samples = 100000U;
    auto const particles = vdf.emit(n_samples);

    auto pos_mom1 = Real{};
    auto pos_mom2 = Real{};
    auto vel_mom1 = Vector{};
    auto vel_mom2 = FourTensor{};
    for (auto const &ptl : particles) {
        pos_mom1 += ptl.pos_x / n_samples;
        pos_mom2 += ptl.pos_x * ptl.pos_x / n_samples;

        vel_mom1 += ptl.vel() / n_samples;

        Scalar const M00 = c * c * ptl.gamma;
        Vector const M0i = c * ptl.g_vel;
        Tensor       Mij;
        Mij.lo() = Mij.hi() = ptl.vel();
        Mij.lo() *= ptl.g_vel;
        Mij.hi() *= { ptl.g_vel.y, ptl.g_vel.z, ptl.g_vel.x };
        vel_mom2 += { M00 / n_samples, M0i / n_samples, Mij / n_samples };
    }
    vel_mom1 = geo.cart2fac(vel_mom1);
    vel_mom2 = geo.cart2fac(vel_mom2);

    auto const X_exact = extent.mean();
    CHECK(std::abs(pos_mom1 - X_exact) < X_exact * 1e-4);
    auto const xx_exact = (extent.min() * extent.min() + extent.min() * extent.max() + extent.max() * extent.max()) / 3.;
    CHECK(std::abs(pos_mom2 - xx_exact) < xx_exact * 1e-4);

    CHECK(vel_mom1.x == Approx{ nV0.x }.epsilon(4e-3));
    CHECK(vel_mom1.y == Approx{ 0 }.margin(1e-3));
    CHECK(vel_mom1.z == Approx{ 0 }.margin(1e-3));

    CHECK(*vel_mom2.tt == Approx{ *nuv0.tt }.epsilon(1e-3));
    CHECK(vel_mom2.ts.x == Approx{ nuv0.ts.x }.epsilon(4e-3));
    CHECK(vel_mom2.ts.y == Approx{ nuv0.ts.y }.margin(4e-3));
    CHECK(vel_mom2.ts.z == Approx{ nuv0.ts.z }.margin(4e-3));
    CHECK(vel_mom2.ss.xx == Approx{ nuv0.ss.xx }.epsilon(6e-3));
    CHECK(vel_mom2.ss.yy == Approx{ nuv0.ss.yy }.epsilon(5e-3));
    CHECK(vel_mom2.ss.zz == Approx{ nuv0.ss.zz }.epsilon(5e-3));
    CHECK(vel_mom2.ss.xy == Approx{ nuv0.ss.xy }.margin(4e-3));
    CHECK(vel_mom2.ss.yz == Approx{ nuv0.ss.yz }.margin(4e-3));
    CHECK(vel_mom2.ss.zx == Approx{ nuv0.ss.zx }.margin(4e-3));

    for (auto const &ptl : particles) {
        CHECK(vdf.weight(ptl) == Approx{ ptl.psd.weight }.epsilon(1e-10));
    }
}

TEST_CASE("Test libPIC::RelativisticVDFVariant::LossconeVDF", "[libPIC::RelativisticVDFVariant::LossconeVDF]")
{
    Real constexpr O0 = 1., c = 4., theta = 30. * M_PI / 180, op = c;
    auto constexpr extent = Range{ 0, 10 } - 1;
    auto constexpr desc   = LossconePlasmaDesc({ { -O0, op }, 10, ShapeOrder::CIC, ParticleScheme::delta_f, .01 },
                                               .1, 2, { .5, .9 }, 1.1);
    auto const geo        = Geometry{ O0, theta };
    auto const vdf        = RelativisticVDFVariant::make(desc, geo, extent, c);

    CHECK(serialize(static_cast<KineticPlasmaDesc const &>(desc)) == serialize(vdf.plasma_desc()));

    auto const n0 = vdf.n0(0);
    CHECK(*n0 == Approx{ 1 }.epsilon(1e-10));

    auto const nV0 = geo.cart2fac(vdf.nV0(0));
    CHECK(nV0.x == Approx{ 1.10000003039345 }.epsilon(1e-7));
    CHECK(nV0.y == Approx{ 0 }.margin(1e-10));
    CHECK(nV0.z == Approx{ 0 }.margin(1e-10));

    auto const nuv0 = geo.cart2fac(vdf.nuv0(0));
    CHECK(nuv0.tt == Approx{ 16.82086014348655 }.epsilon(1e-3));
    CHECK(nuv0.ts.x == Approx{ 4.638779644484074 }.epsilon(1e-4));
    CHECK(nuv0.ts.y == Approx{ 0 }.margin(1e-10));
    CHECK(nuv0.ts.z == Approx{ 0 }.margin(1e-10));
    CHECK(nuv0.ss.xx == Approx{ 1.323094068733219 }.epsilon(1e-4));
    CHECK(nuv0.ss.yy == Approx{ 0.1370451760078678 }.epsilon(1e-3));
    CHECK(nuv0.ss.zz == Approx{ 0.1370451760078678 }.epsilon(1e-3));
    CHECK(nuv0.ss.xy == Approx{ 0 }.margin(1e-10));
    CHECK(nuv0.ss.yz == Approx{ 0 }.margin(1e-10));
    CHECK(nuv0.ss.zx == Approx{ 0 }.margin(1e-10));

    // sampling
    auto const n_samples = 100000U;
    auto const particles = vdf.emit(n_samples);

    auto pos_mom1 = Real{};
    auto pos_mom2 = Real{};
    auto vel_mom1 = Vector{};
    auto vel_mom2 = FourTensor{};
    for (auto const &ptl : particles) {
        pos_mom1 += ptl.pos_x / n_samples;
        pos_mom2 += ptl.pos_x * ptl.pos_x / n_samples;

        vel_mom1 += ptl.vel() / n_samples;

        Scalar const M00 = c * c * ptl.gamma;
        Vector const M0i = c * ptl.g_vel;
        Tensor       Mij;
        Mij.lo() = Mij.hi() = ptl.vel();
        Mij.lo() *= ptl.g_vel;
        Mij.hi() *= { ptl.g_vel.y, ptl.g_vel.z, ptl.g_vel.x };
        vel_mom2 += { M00 / n_samples, M0i / n_samples, Mij / n_samples };
    }
    vel_mom1 = geo.cart2fac(vel_mom1);
    vel_mom2 = geo.cart2fac(vel_mom2);

    auto const X_exact = extent.mean();
    CHECK(std::abs(pos_mom1 - X_exact) < X_exact * 1e-4);
    auto const xx_exact
        = (extent.min() * extent.min() + extent.min() * extent.max() + extent.max() * extent.max())
        / 3.;
    CHECK(std::abs(pos_mom2 - xx_exact) < xx_exact * 1e-4);

    CHECK(vel_mom1.x == Approx{ nV0.x }.epsilon(4e-3));
    CHECK(vel_mom1.y == Approx{ 0 }.margin(1e-3));
    CHECK(vel_mom1.z == Approx{ 0 }.margin(1e-3));

    CHECK(*vel_mom2.tt == Approx{ *nuv0.tt }.epsilon(1e-3));
    CHECK(vel_mom2.ts.x == Approx{ nuv0.ts.x }.epsilon(4e-3));
    CHECK(vel_mom2.ts.y == Approx{ nuv0.ts.y }.margin(4e-3));
    CHECK(vel_mom2.ts.z == Approx{ nuv0.ts.z }.margin(4e-3));
    CHECK(vel_mom2.ss.xx == Approx{ nuv0.ss.xx }.epsilon(7e-3));
    CHECK(vel_mom2.ss.yy == Approx{ nuv0.ss.yy }.epsilon(5e-3));
    CHECK(vel_mom2.ss.zz == Approx{ nuv0.ss.zz }.epsilon(5e-3));
    CHECK(vel_mom2.ss.xy == Approx{ nuv0.ss.xy }.margin(4e-3));
    CHECK(vel_mom2.ss.yz == Approx{ nuv0.ss.yz }.margin(4e-3));
    CHECK(vel_mom2.ss.zx == Approx{ nuv0.ss.zx }.margin(4e-3));

    for (auto const &ptl : particles) {
        CHECK(vdf.weight(ptl) == Approx{ ptl.psd.weight }.epsilon(1e-10));
    }
}
