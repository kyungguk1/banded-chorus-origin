/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "catch2/catch.hpp"

#include <VDF/Particle.h>
#include <cmath>

TEST_CASE("Test common::Particle", "[common::Particle]")
{
    using COMMON_NAMESPACE::Particle;

    Particle ptl;
    CHECK(ptl.vel.fold(true, [](bool lhs, auto rhs) {
        return lhs && std::isnan(rhs);
    }));
    CHECK(std::isnan(ptl.pos_x));
    CHECK(std::isnan(ptl.delta.f));
    CHECK(std::isnan(ptl.delta.w));

    ptl = Particle{ { 1, 2, 3 }, 4 };
    CHECK(ptl.vel.x == 1);
    CHECK(ptl.vel.y == 2);
    CHECK(ptl.vel.z == 3);
    CHECK(ptl.pos_x == 4);
    CHECK(std::isnan(ptl.delta.f));
    CHECK(std::isnan(ptl.delta.w));
}

TEST_CASE("Test common::RelativisticParticle", "[common::RelativisticParticle]")
{
    using Particle = COMMON_NAMESPACE::RelativisticParticle;

    Particle ptl;
    CHECK(ptl.g_vel.fold(true, [](bool lhs, auto rhs) {
        return lhs && std::isnan(rhs);
    }));
    CHECK(std::isnan(ptl.pos_x));
    CHECK(std::isnan(ptl.delta.f));
    CHECK(std::isnan(ptl.delta.w));
    CHECK(std::isnan(ptl.gamma));

    common::Vector v = { 1, 2, 3 };
    double const   c = 5;
    double const   gamma
        = 1 / std::sqrt((1 - std::sqrt(dot(v, v)) / c) * (1 + std::sqrt(dot(v, v) / c)));
    auto const gv = gamma * v;
    ptl           = Particle{ gv, 4, gamma };
    CHECK(ptl.g_vel.x == gv.x);
    CHECK(ptl.g_vel.y == gv.y);
    CHECK(ptl.g_vel.z == gv.z);
    CHECK(ptl.pos_x == 4);
    CHECK(std::isnan(ptl.delta.f));
    CHECK(std::isnan(ptl.delta.w));
    CHECK(ptl.gamma == gamma);
    auto const vel = ptl.vel();
    CHECK(vel.x == v.x);
    CHECK(vel.y == v.y);
    CHECK(vel.z == v.z);
}
