/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "catch2/catch.hpp"

#include <PIC/RelativisticParticle.h>
#include <cmath>

TEST_CASE("Test LibPIC::RelativisticParticle", "[LibPIC::RelativisticParticle]")
{
    using Particle = RelativisticParticle;

    Particle ptl;
    CHECK(ptl.g_vel.fold(true, [](bool lhs, auto rhs) {
        return lhs && std::isnan(rhs);
    }));
    CHECK(std::isnan(ptl.pos.q1));
    CHECK(std::isnan(ptl.psd.weight));
    CHECK(std::isnan(ptl.psd.real_f));
    CHECK(std::isnan(ptl.psd.marker));
    CHECK(std::isnan(ptl.gamma));

    Vector       v = { 1, 2, 3 };
    double const c = 5;
    double const gamma
        = 1 / std::sqrt((1 - std::sqrt(dot(v, v)) / c) * (1 + std::sqrt(dot(v, v) / c)));
    auto const gv = gamma * v;
    ptl           = Particle{ gv, CurviCoord{ 4 }, gamma };
    CHECK(ptl.g_vel.x == gv.x);
    CHECK(ptl.g_vel.y == gv.y);
    CHECK(ptl.g_vel.z == gv.z);
    CHECK(ptl.pos.q1 == 4);
    CHECK(std::isnan(ptl.psd.weight));
    CHECK(std::isnan(ptl.psd.real_f));
    CHECK(std::isnan(ptl.psd.marker));
    CHECK(ptl.gamma == gamma);
    auto const vel = ptl.vel();
    CHECK(vel.x == v.x);
    CHECK(vel.y == v.y);
    CHECK(vel.z == v.z);
}
