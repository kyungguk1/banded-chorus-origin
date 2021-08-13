/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "catch2/catch.hpp"

#include <VDF/LossconeVDF.h>
#include <VDF/MaxwellianVDF.h>

#include <cmath>

TEST_CASE("Test common::MaxwellianVDF", "[common::MaxwellianVDF]")
{
    auto const O0 = 1., theta = 30. * M_PI / 180, op = 4 * O0, c = op;
    auto const extent = Range{ 0, 10 } - 1;
    auto const geo    = Geometry{ O0, theta };
    auto const desc   = BiMaxPlasmaDesc({ { -O0, op }, 10, common::CIC }, .1, 2, -1);
    auto const vdf    = MaxwellianVDF(geo, extent, desc, c);

    CHECK(1 == *vdf.n0(0));
    CHECK(desc.Vd == dot(vdf.nV0(0), geo.e1));
    auto const nvv0 = vdf.nvv0(1) - Tensor{ desc.Vd * desc.Vd, 0, 0, 0, 0, 0 };
    CHECK(std::abs(1 + 2 * desc.T2_T1 - 2 * trace(nvv0) / desc.beta1) < 1e-14);

    auto const n_samples = 100000U;
    auto const particles = vdf.emit(n_samples);

    auto pos_mom1 = Real{};
    auto pos_mom2 = Real{};
    auto vel_mom1 = Vector{};
    auto vel_mom2 = Tensor{};
    for (Particle const &ptl : particles) {
        pos_mom1 += ptl.pos_x / n_samples;
        pos_mom2 += ptl.pos_x * ptl.pos_x / n_samples;
        vel_mom1 += ptl.vel / n_samples;
        Tensor vv;
        vv.lo() = vv.hi() = ptl.vel;
        vv.lo() *= ptl.vel;
        vv.hi() *= { ptl.vel.y, ptl.vel.z, ptl.vel.x };
        vel_mom2 += vv / n_samples;
    }

    auto const X_exact = extent.mean();
    CHECK(std::abs(pos_mom1 - X_exact) < X_exact * 1e-4);
    auto const xx_exact
        = (extent.min() * extent.min() + extent.min() * extent.max() + extent.max() * extent.max())
        / 3.;
    CHECK(std::abs(pos_mom2 - xx_exact) < xx_exact * 1e-4);

    auto const V_exact = vdf.nV0(0);
    auto const diff1   = (V_exact - vel_mom1).fold(Real{}, [](Real const &a, Real const &b) {
        return a + b * b;
      });
    CHECK(std::sqrt(diff1 / 3) < std::abs(desc.Vd) * 1e-3);

    auto const vv_exact = vdf.nvv0(0);
    auto       diff2    = (vv_exact - vel_mom2).fold(Real{}, [](Real const &a, Real const &b) {
        return a + b * b;
             });
    CHECK(std::sqrt(diff2 / 6) < trace(vv_exact) * 1e-3);
}

TEST_CASE("Test common::LossconeVDF::BiMax", "[common::LossconeVDF::BiMax]")
{
    auto const O0 = 1., theta = 30. * M_PI / 180, op = 4 * O0, c = op;
    auto const extent = Range{ 0, 10 } - 1;
    auto const geo    = Geometry{ O0, theta };
    auto const desc   = BiMaxPlasmaDesc({ { -O0, op }, 10, common::CIC }, .1, 2, -1);
    auto const vdf    = LossconeVDF(geo, extent, LossconePlasmaDesc{ desc }, c);

    CHECK(1 == *vdf.n0(0));
    CHECK(desc.Vd == dot(vdf.nV0(0), geo.e1));
    auto const nvv0 = vdf.nvv0(1) - Tensor{ desc.Vd * desc.Vd, 0, 0, 0, 0, 0 };
    CHECK(std::abs(1 + 2 * desc.T2_T1 - 2 * trace(nvv0) / desc.beta1) < 1e-14);

    auto const n_samples = 100000U;
    auto const particles = vdf.emit(n_samples);

    auto pos_mom1 = Real{};
    auto pos_mom2 = Real{};
    auto vel_mom1 = Vector{};
    auto vel_mom2 = Tensor{};
    for (Particle const &ptl : particles) {
        pos_mom1 += ptl.pos_x / n_samples;
        pos_mom2 += ptl.pos_x * ptl.pos_x / n_samples;
        vel_mom1 += ptl.vel / n_samples;
        Tensor vv;
        vv.lo() = vv.hi() = ptl.vel;
        vv.lo() *= ptl.vel;
        vv.hi() *= { ptl.vel.y, ptl.vel.z, ptl.vel.x };
        vel_mom2 += vv / n_samples;
    }

    auto const X_exact = extent.mean();
    CHECK(std::abs(pos_mom1 - X_exact) < X_exact * 1e-4);
    auto const xx_exact
        = (extent.min() * extent.min() + extent.min() * extent.max() + extent.max() * extent.max())
        / 3.;
    CHECK(std::abs(pos_mom2 - xx_exact) < xx_exact * 1e-4);

    auto const V_exact = vdf.nV0(0);
    auto const diff1   = (V_exact - vel_mom1).fold(Real{}, [](Real const &a, Real const &b) {
        return a + b * b;
      });
    CHECK(std::sqrt(diff1 / 3) < std::abs(desc.Vd) * 1e-3);

    auto const vv_exact = vdf.nvv0(0);
    auto       diff2    = (vv_exact - vel_mom2).fold(Real{}, [](Real const &a, Real const &b) {
        return a + b * b;
             });
    CHECK(std::sqrt(diff2 / 6) < trace(vv_exact) * 1e-3);
}

TEST_CASE("Test common::LossconeVDF::Loss", "[common::LossconeVDF::Loss]")
{
    auto const O0 = 1., theta = 30. * M_PI / 180, op = 4 * O0, c = op;
    auto const extent = Range{ 0, 10 } - 1;
    auto const geo    = Geometry{ O0, theta };
    auto const desc = LossconePlasmaDesc({ { -O0, op }, 10, common::CIC }, .1, 2, { .5, .9 }, 1.1);
    auto const vdf  = LossconeVDF(geo, extent, desc, c);

    CHECK(1 == *vdf.n0(0));
    CHECK(desc.Vd == dot(vdf.nV0(0), geo.e1));
    auto const nvv0 = vdf.nvv0(1) - Tensor{ desc.Vd * desc.Vd, 0, 0, 0, 0, 0 };
    CHECK(std::abs(1 + 2 * desc.T2_T1 - 2 * trace(nvv0) / desc.beta1) < 1e-14);

    auto const n_samples = 100000U;
    auto const particles = vdf.emit(n_samples);

    auto pos_mom1 = Real{};
    auto pos_mom2 = Real{};
    auto vel_mom1 = Vector{};
    auto vel_mom2 = Tensor{};
    for (Particle const &ptl : particles) {
        pos_mom1 += ptl.pos_x / n_samples;
        pos_mom2 += ptl.pos_x * ptl.pos_x / n_samples;
        vel_mom1 += ptl.vel / n_samples;
        Tensor vv;
        vv.lo() = vv.hi() = ptl.vel;
        vv.lo() *= ptl.vel;
        vv.hi() *= { ptl.vel.y, ptl.vel.z, ptl.vel.x };
        vel_mom2 += vv / n_samples;
    }

    auto const X_exact = extent.mean();
    CHECK(std::abs(pos_mom1 - X_exact) < X_exact * 1e-4);
    auto const xx_exact
        = (extent.min() * extent.min() + extent.min() * extent.max() + extent.max() * extent.max())
        / 3.;
    CHECK(std::abs(pos_mom2 - xx_exact) < xx_exact * 1e-4);

    auto const V_exact = vdf.nV0(0);
    auto const diff1   = (V_exact - vel_mom1).fold(Real{}, [](Real const &a, Real const &b) {
        return a + b * b;
      });
    CHECK(std::sqrt(diff1 / 3) < std::abs(desc.Vd) * 1e-3);

    auto const vv_exact = vdf.nvv0(0);
    auto       diff2    = (vv_exact - vel_mom2).fold(Real{}, [](Real const &a, Real const &b) {
        return a + b * b;
             });
    CHECK(std::sqrt(diff2 / 6) < trace(vv_exact) * 1e-3);
}
