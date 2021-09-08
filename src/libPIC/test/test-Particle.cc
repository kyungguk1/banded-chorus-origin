/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "catch2/catch.hpp"

#include <PIC/Particle.h>
#include <cmath>

TEST_CASE("Test libPIC::Particle", "[libPIC::Particle]")
{
    Particle ptl;
    CHECK(ptl.vel.fold(true, [](bool lhs, auto rhs) {
        return lhs && std::isnan(rhs);
    }));
    CHECK(std::isnan(ptl.pos_x));
    CHECK(std::isnan(ptl.psd.f));
    CHECK(std::isnan(ptl.psd.w));
    CHECK(-1 == ptl.id);

    for (long i = 0; i < 100; ++i) {
        ptl = Particle{ { 1, 2, 3 }, 4 };
        (void)ptl;
    }
    ptl = Particle{ { 1, 2, 3 }, 4 };
    CHECK(ptl.vel.x == 1);
    CHECK(ptl.vel.y == 2);
    CHECK(ptl.vel.z == 3);
    CHECK(ptl.pos_x == 4);
    CHECK(std::isnan(ptl.psd.f));
    CHECK(std::isnan(ptl.psd.w));
    // CHECK(100 == ptl.id);
}
