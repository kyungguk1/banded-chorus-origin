/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "catch2/catch.hpp"

#include <PIC/MirrorGeometry.h>

TEST_CASE("Test libPIC::MirrorGeometry", "[libPIC::MirrorGeometry]")
{
    CHECK_THROWS_AS(MirrorGeometry(-1, 1), std::invalid_argument);
    CHECK_THROWS_AS(MirrorGeometry(1, 0), std::invalid_argument);
    CHECK_NOTHROW(MirrorGeometry(0, 1));

    { // homogeneous
        constexpr Real       xi = 0;
        constexpr Real       D1 = 0.1;
        MirrorGeometry const mirror{ xi, D1 };

        CHECK(xi == mirror.xi());
        CHECK(xi * xi == mirror.xi2());
        CHECK(mirror.is_homogeneous());

        CHECK(mirror.D().x == D1);
        CHECK(mirror.D().y == 1);
        CHECK(mirror.D().z == 1);
        CHECK(mirror.D1() == D1);
        CHECK(mirror.D2() == 1);
        CHECK(mirror.D3() == 1);

        CHECK(mirror.sqrt_g() == D1);
        CHECK(mirror.det_gij() == D1 * D1);

        CHECK(mirror.is_valid(CurviCoord{ 0 }));
        CHECK(mirror.is_valid(CurviCoord{ 1 }));
        CHECK(mirror.is_valid(CurviCoord{ -100 }));
    }

    { // inhomogeneous
        constexpr Real       xi = 0.112;
        constexpr Real       D1 = 2;
        MirrorGeometry const mirror{ xi, D1 };

        CHECK(xi == mirror.xi());
        CHECK(xi * xi == mirror.xi2());
        CHECK(!mirror.is_homogeneous());

        CHECK(mirror.D().x == D1);
        CHECK(mirror.D().y == 1);
        CHECK(mirror.D().z == 1);
        CHECK(mirror.D1() == D1);
        CHECK(mirror.D2() == 1);
        CHECK(mirror.D3() == 1);

        CHECK(mirror.sqrt_g() == D1);
        CHECK(mirror.det_gij() == D1 * D1);

        CHECK(mirror.is_valid(CurviCoord{ 0 }));
        CHECK(mirror.is_valid(CurviCoord{ M_PI_2 * 0.99999999 / (xi * D1) }));
        CHECK(!mirror.is_valid(CurviCoord{ -M_PI_2 * 1.00000001 / (xi * D1) }));
    }
}

TEST_CASE("Test libPIC::MirrorGeometry::Cotrans", "[libPIC::MirrorGeometry::Cotrans]")
{
    { // homogeneous
        constexpr Real       xi = 0;
        constexpr Real       D1 = 0.1;
        MirrorGeometry const mirror{ xi, D1 };

        constexpr CartCoord cart1{ 1.45 };
        auto const          curvi = mirror.cotrans(cart1);
        CHECK(curvi.q1 * D1 == Approx{ cart1.x }.epsilon(1e-10));
        auto const cart2 = mirror.cotrans(curvi);
        CHECK(cart2.x == Approx{ cart1.x }.epsilon(1e-10));
    }

    { // inhomogeneous
        constexpr Real       xi = 0.112;
        constexpr Real       D1 = 2;
        MirrorGeometry const mirror{ xi, D1 };

        constexpr CartCoord cart1{ 1.45 };
        auto const          curvi = mirror.cotrans(cart1);
        CHECK(curvi.q1 * D1 == Approx{ 1.4374506757616818 }.epsilon(1e-10));
        auto const cart2 = mirror.cotrans(curvi);
        CHECK(cart2.x == Approx{ cart1.x }.epsilon(1e-10));
    }
}
