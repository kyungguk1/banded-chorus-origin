/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "catch2/catch.hpp"

#include <VDF/Particle.h>
#include <cmath>

using COMMON_NAMESPACE::Particle;

TEST_CASE("Test common::Particle", "[common::Particle]")
{
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
