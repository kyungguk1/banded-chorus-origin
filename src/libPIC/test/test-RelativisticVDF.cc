/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "catch2/catch.hpp"

#include <PIC/RelativisticMaxwellianVDF.h>
#include <PIC/RelativisticVDF.h>

#include <array>
#include <cmath>
#include <vector>

TEST_CASE("Test libPIC::RelativisticMaxwellianVDF", "[libPIC::RelativisticMaxwellianVDF]")
{
    Real constexpr O0 = 1., c = 4., theta = 30. * M_PI / 180, op = c;
    auto constexpr extent = Range{ 0, 10 } - 1;
    auto constexpr desc   = BiMaxPlasmaDesc(
          { { -O0, op }, 10, ShapeOrder::CIC, ParticleScheme::delta_f, .01 }, .1, 2, -1);
    auto const geo = Geometry{ O0, theta };
    auto const vdf = RelativisticMaxwellianVDF(desc, geo, extent, c);

    CHECK(serialize(desc) == serialize(vdf.plasma_desc()));

    // check equilibrium distribution
    std::vector<std::array<Real, 4>> const psd_fac{
        { -1.6392524206218977, 0.021875942164225923, 0.04925485127062607, 0.0973188880413619 },
        { -1.5051260579112105, 0.040148281187166096, 0.05150777282466861, 0.35247699708249824 },
        { -0.7830701882394062, 0.005700257970058733, -0.1172865493600963, 1.414314963192856 },
        { -0.9785033720943546, -0.49923213712571507, -0.10370922507536744, 0.7221584063296483 },
        { -1.350221774077486, 0.2713083165011999, -0.09889433254807288, 0.7282911713582314 },
        { -0.9508297875190594, 0.37508038650353426, -0.313592326586365, 0.7714216315891262 },
        { -0.8336099327629195, 1.1659935984594576, -0.4128767284001784, 0.0007310532694093008 },
        { -1.0564823801259944, 0.23153027092986664, 1.1613413752164727, 0.002468299778026409 },
        { -0.9017292106907236, -0.05702224768678881, 1.2069909926527644, 0.0013790903785312226 },
        { -0.8436544349331425, 0.48465237345520956, 0.9111339599881719, 0.0083446540722853 },
        { -0.8652798796112622, 0.11048137794206545, 0.09882412895049376, 1.884125958696743 },
        { -1.7554071131177433, -0.450422767876925, 0.9683991879664678, 0.0001308425167868628 },
        { -1.9428271847011391, 0.2738650192628952, -0.11212528492650621, 0.0011697835488729265 },
        { -1.972413772006004, -0.014908947632139463, 1.281505591428828, 6.543249366634526e-7 },
        { -0.4818929955170441, -0.3927796297731702, -1.1932628317966936, 0.0000309359260620754 },
        { -1.6765890166556645, 0.7939876808249722, -0.369216657432603, 0.0018491046538197107 },
        { -1.5486701554929136, 0.049592709815474596, 0.6571856387149972, 0.031605006791340394 },
        { -1.5756670313727115, -0.0859380045884064, -1.0909884545727935, 0.0006761005635692306 },
        { -1.2002311116799431, -0.36096909371130254, 0.3377820771662346, 0.6394894974185505 },
        { -1.8008368913564563, -0.13369676439224762, -0.7986313984991514, 0.0006881792319621204 }
    };
    for (auto const &vals : psd_fac) {
        auto const &[u1, u2, u3, psd1] = vals;
        auto const ptl                 = RelativisticParticle(geo.fac2cart({ u1, u2, u3 }), 0, 0);
        auto const psd2                = vdf.f0(ptl);
        auto const psd3                = vdf.g0(ptl);

        INFO("u1 = " << u1 << ", u2 = " << u2 << ", u3 = " << u3 << ", psd1 = " << psd1
                     << ", psd2 = " << psd2 << ", psd3 = " << psd3);
        CHECK(psd2 == Approx{ psd1 }.epsilon(1e-10));
        CHECK(psd3 == Approx{ psd1 }.epsilon(1e-10));
    }

    // check equilibrium macro variables
    CHECK(vdf.n(0) / vdf.n0(0)
          == Approx{ c / std::sqrt((c - desc.Vd) * (c + desc.Vd)) }.epsilon(1e-10));

    auto const n00 = *vdf.n00c2(0) / (c * c);
    CHECK(n00 == Approx{ 0.975760381496556 }.epsilon(1e-4));
    CHECK(n00 == Approx{ 0.97581025402636 }.epsilon(1e-7));

    auto const P0Om0 = vdf.P0Om0(0);
    CHECK(P0Om0.xx == Approx{ 0.04789494769856955 }.epsilon(1e-3));
    CHECK(P0Om0.xx == Approx{ 0.04788278246309458 }.epsilon(1e-7));

    auto const nU0 = geo.cart2fac(vdf.nU0(0));
    CHECK(nU0.t == Approx{ 4.131182252881983 }.epsilon(1e-8));
    CHECK(nU0.s.x == Approx{ -1.032795557455364 }.epsilon(1e-8));
    CHECK(nU0.s.y == Approx{ 0 }.margin(1e-10));
    CHECK(nU0.s.z == Approx{ 0 }.margin(1e-10));

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
    auto vel_mom1 = FourVector{};
    auto vel_mom2 = FourTensor{};
    for (RelativisticParticle const &ptl : particles) {
        pos_mom1 += ptl.pos_x / n_samples;
        pos_mom2 += ptl.pos_x * ptl.pos_x / n_samples;

        vel_mom1 += { c / n_samples, ptl.vel() / n_samples };

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

    auto const nV0 = nU0 * std::sqrt((1 - desc.Vd / c) * (1 + desc.Vd / c));
    CHECK(*vel_mom1.t == Approx{ *nV0.t }.epsilon(1e-10));
    CHECK(vel_mom1.s.x == Approx{ nV0.s.x }.epsilon(4e-3));
    CHECK(vel_mom1.s.y == Approx{ 0 }.margin(1e-3));
    CHECK(vel_mom1.s.z == Approx{ 0 }.margin(1e-3));

    CHECK(*vel_mom2.tt == Approx{ *nuv0.tt }.epsilon(1e-3));
    CHECK(vel_mom2.ts.x == Approx{ nuv0.ts.x }.epsilon(4e-3));
    CHECK(vel_mom2.ts.y == Approx{ nuv0.ts.y }.margin(4e-3));
    CHECK(vel_mom2.ts.z == Approx{ nuv0.ts.z }.margin(4e-3));
    CHECK(vel_mom2.ss.xx == Approx{ nuv0.ss.xx }.epsilon(5e-3));
    CHECK(vel_mom2.ss.yy == Approx{ nuv0.ss.yy }.epsilon(5e-3));
    CHECK(vel_mom2.ss.zz == Approx{ nuv0.ss.zz }.epsilon(5e-3));
    CHECK(vel_mom2.ss.xy == Approx{ nuv0.ss.xy }.margin(4e-3));
    CHECK(vel_mom2.ss.yz == Approx{ nuv0.ss.yz }.margin(4e-3));
    CHECK(vel_mom2.ss.zx == Approx{ nuv0.ss.zx }.margin(4e-3));

    for (RelativisticParticle const &ptl : particles) {
        REQUIRE(vdf.weight(ptl) == Approx{ ptl.psd.w }.epsilon(1e-10));
    }
}
