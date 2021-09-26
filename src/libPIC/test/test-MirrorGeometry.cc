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

        CHECK(mirror.inv_D().x == 1 / D1);
        CHECK(mirror.inv_D().y == 1 / 1);
        CHECK(mirror.inv_D().z == 1 / 1);
        CHECK(mirror.inv_D1() == 1 / D1);
        CHECK(mirror.inv_D2() == 1 / 1);
        CHECK(mirror.inv_D3() == 1 / 1);

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

        constexpr CartCoord cart1{ 14.5 };
        auto const          curvi = mirror.cotrans(cart1);
        CHECK(curvi.q1 * D1 == Approx{ cart1.x }.epsilon(1e-10));
        auto const cart2 = mirror.cotrans(curvi);
        CHECK(cart2.x == Approx{ cart1.x }.epsilon(1e-10));
    }

    { // inhomogeneous
        constexpr Real       xi = 0.112;
        constexpr Real       D1 = 2;
        MirrorGeometry const mirror{ xi, D1 };

        constexpr CartCoord cart1{ 14.5 };
        auto const          curvi = mirror.cotrans(cart1);
        CHECK(curvi.q1 * D1 == Approx{ 9.09702270985558 }.epsilon(1e-10));
        auto const cart2 = mirror.cotrans(curvi);
        CHECK(cart2.x == Approx{ cart1.x }.epsilon(1e-10));
    }
}

TEST_CASE("Test libPIC::MirrorGeometry::Field", "[libPIC::MirrorGeometry::Field]")
{
    { // homogeneous
        constexpr Real       xi = 0;
        constexpr Real       D1 = 0.1;
        MirrorGeometry const mirror{ xi, D1 };

        constexpr CartCoord cart{ 14.5 };
        auto const          curvi = mirror.cotrans(cart);
        Vector              B;

        B = mirror.Bcart_div_B0(cart);
        CHECK(B.x == 1);
        CHECK(B.y == 0);
        CHECK(B.z == 0);

        B = mirror.Bcart_div_B0(curvi);
        CHECK(B.x == 1);
        CHECK(B.y == 0);
        CHECK(B.z == 0);

        B = mirror.Bcart_div_B0(cart, 10, 20);
        CHECK(B.x == 1);
        CHECK(B.y == 0);
        CHECK(B.z == 0);

        B = mirror.Bcart_div_B0(curvi, 10, 20);
        CHECK(B.x == 1);
        CHECK(B.y == 0);
        CHECK(B.z == 0);

        B = mirror.Bcontr_div_B0(cart);
        CHECK(B.x * D1 == Approx{ 1 }.epsilon(1e-10));
        CHECK(B.y == 0);
        CHECK(B.z == 0);

        B = mirror.Bcontr_div_B0(curvi);
        CHECK(B.x * D1 == Approx{ 1 }.epsilon(1e-10));
        CHECK(B.y == 0);
        CHECK(B.z == 0);

        B = mirror.Bcovar_div_B0(cart);
        CHECK(B.x / D1 == Approx{ 1 }.epsilon(1e-10));
        CHECK(B.y == 0);
        CHECK(B.z == 0);

        B = mirror.Bcovar_div_B0(curvi);
        CHECK(B.x / D1 == Approx{ 1 }.epsilon(1e-10));
        CHECK(B.y == 0);
        CHECK(B.z == 0);
    }

    { // inhomogeneous
        constexpr Real       xi = 0.112;
        constexpr Real       D1 = 2;
        MirrorGeometry const mirror{ xi, D1 };

        constexpr CartCoord cart{ 14.5 };
        auto const          curvi = mirror.cotrans(cart);
        Vector              B;

        B = mirror.Bcart_div_B0(cart);
        CHECK(B.x == Approx{ 3.637376 }.epsilon(1e-10));
        CHECK(B.y == 0);
        CHECK(B.z == 0);

        B = mirror.Bcart_div_B0(curvi);
        CHECK(B.x == Approx{ 3.637376 }.epsilon(1e-10));
        CHECK(B.y == 0);
        CHECK(B.z == 0);

        B = mirror.Bcart_div_B0(cart, 10, 20);
        CHECK(B.x == Approx{ 3.637376 }.epsilon(1e-10));
        CHECK(B.y == Approx{ -1.8188800000000003 }.epsilon(1e-10));
        CHECK(B.z == Approx{ -3.6377600000000005 }.epsilon(1e-10));

        B = mirror.Bcart_div_B0(curvi, 10, 20);
        CHECK(B.x == Approx{ 3.637376 }.epsilon(1e-10));
        CHECK(B.y == Approx{ -1.8188800000000003 }.epsilon(1e-10));
        CHECK(B.z == Approx{ -3.6377600000000005 }.epsilon(1e-10));

        B = mirror.Bcontr_div_B0(cart);
        CHECK(B.x * D1 == Approx{ 1 }.epsilon(1e-10));
        CHECK(B.y == 0);
        CHECK(B.z == 0);

        B = mirror.Bcontr_div_B0(curvi);
        CHECK(B.x * D1 == Approx{ 1 }.epsilon(1e-10));
        CHECK(B.y == 0);
        CHECK(B.z == 0);

        B = mirror.Bcovar_div_B0(cart);
        CHECK(B.x / D1 == Approx{ dot(mirror.Bcart_div_B0(cart), mirror.Bcart_div_B0(cart)) }.epsilon(1e-10));
        CHECK(B.y == 0);
        CHECK(B.z == 0);

        B = mirror.Bcovar_div_B0(curvi);
        CHECK(B.x / D1 == Approx{ dot(mirror.Bcart_div_B0(cart), mirror.Bcart_div_B0(cart)) }.epsilon(1e-10));
        CHECK(B.y == 0);
        CHECK(B.z == 0);
    }
}
