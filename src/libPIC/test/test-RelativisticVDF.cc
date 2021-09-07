/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "catch2/catch.hpp"

#include <PIC/RelativisticLossconeVDF.h>
#include <PIC/RelativisticMaxwellianVDF.h>
#include <PIC/RelativisticVDF.h>

#include <array>
#include <cmath>
#include <vector>

TEST_CASE("Test libPIC::RelativisticMaxwellianVDF", "[libPIC::RelativisticMaxwellianVDF]")
{
    Real constexpr O0 = 1., c = 4., theta = 30. * M_PI / 180, op = c;
    auto constexpr extent = Range{ 0, 10 } - 1;
    auto constexpr desc   = BiMaxPlasmaDesc({ { -O0, op }, 10, ShapeOrder::CIC, ParticleScheme::delta_f, .01 },
                                            .1, 2, -1);
    auto const geo        = Geometry{ O0, theta };
    auto const vdf        = RelativisticMaxwellianVDF(desc, geo, extent, c);

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

        INFO("u1 = " << u1 << ", u2 = " << u2 << ", u3 = " << u3 << ", psd1 = " << psd1 << ", psd2 = " << psd2 << ", psd3 = " << psd3);
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
        CHECK(vdf.weight(ptl) == Approx{ ptl.psd.w }.epsilon(1e-10));
    }
}

TEST_CASE("Test libPIC::RelativisticLossconeVDF::BiMax", "[libPIC::RelativisticLossconeVDF::BiMax]")
{
    Real constexpr O0 = 1., c = 4., theta = 30. * M_PI / 180, op = c;
    auto constexpr extent = Range{ 0, 10 } - 1;
    auto constexpr desc   = BiMaxPlasmaDesc({ { -O0, op }, 10, ShapeOrder::CIC, ParticleScheme::delta_f, .01 },
                                            .1, 2, -1);
    auto const geo        = Geometry{ O0, theta };
    auto const vdf        = RelativisticLossconeVDF(LossconePlasmaDesc{ desc }, geo, extent, c);

    CHECK(serialize(LossconePlasmaDesc{ desc }) == serialize(vdf.plasma_desc()));

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

        INFO("u1 = " << u1 << ", u2 = " << u2 << ", u3 = " << u3 << ", psd1 = " << psd1 << ", psd2 = " << psd2 << ", psd3 = " << psd3);
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
        CHECK(vdf.weight(ptl) == Approx{ ptl.psd.w }.epsilon(1e-10));
    }
}

TEST_CASE("Test libPIC::RelativisticLossconeVDF::Losscone", "[libPIC::RelativisticLossconeVDF::Losscone]")
{
    Real constexpr O0 = 1., c = 4., theta = 30. * M_PI / 180, op = c;
    auto constexpr extent = Range{ 0, 10 } - 1;
    auto constexpr desc   = LossconePlasmaDesc({ { -O0, op }, 10, ShapeOrder::CIC, ParticleScheme::delta_f, .01 },
                                               .1, 2, { .5, .9 }, 1.1);
    auto const geo        = Geometry{ O0, theta };
    auto const vdf        = RelativisticLossconeVDF(desc, geo, extent, c);

    CHECK(serialize(desc) == serialize(vdf.plasma_desc()));

    // check equilibrium distribution
    std::vector<std::array<Real, 4>> const psd_fac{
        { 1.830946816273993, 1.101015375601375, -0.02782941713257447, 0.0004684370166405492 },
        { 0.2457196872446087, 0.03695465231273989, 0.2290062025646392, 0.0004678623456198146 },
        { 1.906469500981201, -0.1168718362495395, 0.6966531923004128, 0.0028890947308751 },
        { 1.410843230405879, 0.918983440600655, -0.903737919893139, 0.001591125764730341 },
        { 2.060879021176165, 0.1517422583845632, -0.03718902345982182, 0.000905588123638418 },
        { 2.1641464916788, 0.00542836536090446, -0.2754324192784612, 0.0001658609885420265 },
        { 1.103354949008349, -1.088211371270094, 0.6557962996982559, 0.002703385900777482 },
        { 1.229805957952653, -0.1473395909430888, 0.2065024831332222, 1.254155588046163 },
        { 1.855271880403062, 0.6382355211857722, 0.3270508751861647, 0.00518644908888417 },
        { 1.985436851502545, -0.0823254187172925, -0.06943195229628792, 0.002782661981653955 },
        { 1.592811295429149, -0.7048792906653926, -0.4240924036554178, 0.03775237561821107 },
        { 0.890457042928043, -0.5064793127784148, -0.1012665072553588, 0.4457174481983278 },
        { 0.951048152198153, -0.04013633404371774, 0.803771554604695, 0.1390742901845188 },
        { 0.4840716107513686, 0.4517697851347681, -0.3487992159234443, 0.00906483373247185 },
        { 0.6502174217257453, -1.144991740714791, 0.3859752433325411, 0.0003519520565018298 },
        { 0.3495978335783963, 0.6054963677309365, 0.6284613341925716, 0.0001855909741311024 },
        { 1.782126610537641, -1.277017482718514, -0.2908740770393852, 0.0000929673966632176 },
        { 0.993612266179765, -0.02834690460920546, 0.00916011794555343, 1.105257797864694 },
        { 1.227303231797827, -0.3147884765773692, -0.3338675691216383, 0.950750174451146 },
        { 0.910299633196709, 0.2525480653188529, -0.7443207674158588, 0.1313271901697628 }
    };
    for (auto const &vals : psd_fac) {
        auto const &[u1, u2, u3, psd1] = vals;
        auto const ptl                 = RelativisticParticle(geo.fac2cart({ u1, u2, u3 }), 0, 0);
        auto const psd2                = vdf.f0(ptl);
        auto const psd3                = vdf.g0(ptl);

        INFO("u1 = " << u1 << ", u2 = " << u2 << ", u3 = " << u3 << ", psd1 = " << psd1 << ", psd2 = " << psd2 << ", psd3 = " << psd3);
        CHECK(psd2 == Approx{ psd1 }.epsilon(1e-10));
        CHECK(psd3 == Approx{ psd1 }.epsilon(1e-10));
    }

    // check equilibrium macro variables
    CHECK(vdf.n(0) / vdf.n0(0)
          == Approx{ c / std::sqrt((c - desc.Vd) * (c + desc.Vd)) }.epsilon(1e-10));

    auto const n00 = *vdf.n00c2(0) / (c * c);
    CHECK(n00 == Approx{ 0.971574683500847 }.epsilon(1e-4));
    CHECK(n00 == Approx{ 0.971659543310232 }.epsilon(1e-7));

    auto const P0Om0 = vdf.P0Om0(0);
    CHECK(P0Om0.xx == Approx{ 0.04742967300474676 }.epsilon(1e-3));
    CHECK(P0Om0.xx == Approx{ 0.04741121858556467 }.epsilon(1e-7));

    auto const nU0 = geo.cart2fac(vdf.nU0(0));
    CHECK(nU0.t == Approx{ 4.160407802694631 }.epsilon(1e-7));
    CHECK(nU0.s.x == Approx{ 1.144112140679263 }.epsilon(1e-7));
    CHECK(nU0.s.y == Approx{ 0 }.margin(1e-10));
    CHECK(nU0.s.z == Approx{ 0 }.margin(1e-10));

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
    CHECK(vel_mom2.ss.xx == Approx{ nuv0.ss.xx }.epsilon(7e-3));
    CHECK(vel_mom2.ss.yy == Approx{ nuv0.ss.yy }.epsilon(5e-3));
    CHECK(vel_mom2.ss.zz == Approx{ nuv0.ss.zz }.epsilon(5e-3));
    CHECK(vel_mom2.ss.xy == Approx{ nuv0.ss.xy }.margin(4e-3));
    CHECK(vel_mom2.ss.yz == Approx{ nuv0.ss.yz }.margin(4e-3));
    CHECK(vel_mom2.ss.zx == Approx{ nuv0.ss.zx }.margin(4e-3));

    for (RelativisticParticle const &ptl : particles) {
        CHECK(vdf.weight(ptl) == Approx{ ptl.psd.w }.epsilon(1e-10));
    }
}
