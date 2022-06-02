/*
 * Copyright (c) 2021-2022, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "catch2/catch.hpp"

#define LIBPIC_INLINE_VERSION 1
#include <PIC/Geometry/MirrorGeometry.h>

using Detail::MirrorGeometry;

TEST_CASE("Test LibPIC::MirrorGeometry", "[LibPIC::MirrorGeometry]")
{
    CHECK_THROWS_AS(MirrorGeometry(-1, 1), std::invalid_argument);
    CHECK_THROWS_AS(MirrorGeometry(1, 0), std::invalid_argument);
    CHECK_NOTHROW(MirrorGeometry(0, 1));

    { // homogeneous
        constexpr Real       xi = 0;
        constexpr Real       D1 = 0.1;
        constexpr Real       D2 = 0.43;
        constexpr Real       D3 = 1.54;
        MirrorGeometry const mirror{ xi, { D1, D2, D3 } };

        CHECK(xi == mirror.xi());
        CHECK(xi * xi == mirror.xi2());
        CHECK(mirror.is_homogeneous());

        CHECK(mirror.D().x == D1);
        CHECK(mirror.D().y == D2);
        CHECK(mirror.D().z == D3);
        CHECK(mirror.D1() == D1);
        CHECK(mirror.D2() == D2);
        CHECK(mirror.D3() == D3);

        CHECK(mirror.inv_D().x == 1 / D1);
        CHECK(mirror.inv_D().y == 1 / D2);
        CHECK(mirror.inv_D().z == 1 / D3);
        CHECK(mirror.inv_D1() == 1 / D1);
        CHECK(mirror.inv_D2() == 1 / D2);
        CHECK(mirror.inv_D3() == 1 / D3);

        CHECK(mirror.sqrt_g() == D1 * D2 * D3);
        CHECK(mirror.det_gij() == mirror.sqrt_g() * mirror.sqrt_g());

        CHECK(mirror.is_valid(CurviCoord{ 0 }));
        CHECK(mirror.is_valid(CurviCoord{ 1 }));
        CHECK(mirror.is_valid(CurviCoord{ -100 }));
    }

    { // inhomogeneous
        constexpr Real       xi = 0.112;
        constexpr Real       D1 = 2;
        constexpr Real       D2 = 0.43;
        constexpr Real       D3 = 1.54;
        MirrorGeometry const mirror{ xi, { D1, D2, D3 } };

        CHECK(xi == mirror.xi());
        CHECK(xi * xi == mirror.xi2());
        CHECK(!mirror.is_homogeneous());

        CHECK(mirror.D().x == D1);
        CHECK(mirror.D().y == D2);
        CHECK(mirror.D().z == D3);
        CHECK(mirror.D1() == D1);
        CHECK(mirror.D2() == D2);
        CHECK(mirror.D3() == D3);

        CHECK(mirror.inv_D().x == 1 / D1);
        CHECK(mirror.inv_D().y == 1 / D2);
        CHECK(mirror.inv_D().z == 1 / D3);
        CHECK(mirror.inv_D1() == 1 / D1);
        CHECK(mirror.inv_D2() == 1 / D2);
        CHECK(mirror.inv_D3() == 1 / D3);

        CHECK(mirror.sqrt_g() == D1 * D2 * D3);
        CHECK(mirror.det_gij() == mirror.sqrt_g() * mirror.sqrt_g());

        CHECK(mirror.is_valid(CurviCoord{ 0 }));
        CHECK(mirror.is_valid(CurviCoord{ M_PI_2 * 0.99999999 / (xi * D1) }));
        CHECK(!mirror.is_valid(CurviCoord{ -M_PI_2 * 1.00000001 / (xi * D1) }));
    }
}

TEST_CASE("Test LibPIC::MirrorGeometry::Cotrans", "[LibPIC::MirrorGeometry::Cotrans]")
{
    { // homogeneous
        constexpr Real       xi = 0;
        constexpr Real       D1 = 0.1;
        constexpr Real       D2 = 0.43;
        constexpr Real       D3 = 1.54;
        MirrorGeometry const mirror{ xi, { D1, D2, D3 } };

        constexpr CartCoord cart1{ 14.5 };
        auto const          curvi = mirror.cotrans(cart1);
        CHECK(curvi.q1 * D1 == Approx{ cart1.x }.epsilon(1e-10));
        auto const cart2 = mirror.cotrans(curvi);
        CHECK(cart2.x == Approx{ cart1.x }.epsilon(1e-10));
    }

    { // inhomogeneous
        constexpr Real       xi = 0.112;
        constexpr Real       D1 = 2;
        constexpr Real       D2 = 0.43;
        constexpr Real       D3 = 1.54;
        MirrorGeometry const mirror{ xi, { D1, D2, D3 } };

        constexpr CartCoord cart1{ 14.5 };
        auto const          curvi = mirror.cotrans(cart1);
        CHECK(curvi.q1 * D1 == Approx{ 9.09702270985558 }.epsilon(1e-10));
        auto const cart2 = mirror.cotrans(curvi);
        CHECK(cart2.x == Approx{ cart1.x }.epsilon(1e-10));
    }
}

TEST_CASE("Test LibPIC::MirrorGeometry::Field", "[LibPIC::MirrorGeometry::Field]")
{
    { // homogeneous
        constexpr Real       xi = 0;
        constexpr Real       D1 = 0.1;
        constexpr Real       D2 = 0.43;
        constexpr Real       D3 = 1.54;
        MirrorGeometry const mirror{ xi, { D1, D2, D3 } };

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

        CHECK(std::pow(mirror.Bmag_div_B0(cart), 2) == Approx{ dot(B, B) }.epsilon(1e-10));
        CHECK(std::pow(mirror.Bmag_div_B0(curvi), 2) == Approx{ dot(B, B) }.epsilon(1e-10));

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
        constexpr Real       D2 = 0.43;
        constexpr Real       D3 = 1.54;
        MirrorGeometry const mirror{ xi, { D1, D2, D3 } };

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

        CHECK(std::pow(mirror.Bmag_div_B0(cart), 2) == Approx{ dot(B, B) }.epsilon(1e-10));
        CHECK(std::pow(mirror.Bmag_div_B0(curvi), 2) == Approx{ dot(B, B) }.epsilon(1e-10));

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

TEST_CASE("Test LibPIC::MirrorGeometry::Basis", "[LibPIC::MirrorGeometry::Basis]")
{
    { // homogeneous
        constexpr Real       xi = 0;
        constexpr Real       D1 = 0.1;
        constexpr Real       D2 = 0.43;
        constexpr Real       D3 = 1.54;
        MirrorGeometry const mirror{ xi, { D1, D2, D3 } };

        constexpr CartCoord cart{ 14.5 };
        auto const          curvi = mirror.cotrans(cart);
        Vector              basis;
        Tensor              bases;

        // covar metric
        bases = mirror.covar_metric(cart);
        CHECK(bases.xx == Approx{ D1 * D1 }.epsilon(1e-10));
        CHECK(bases.yy == Approx{ D2 * D2 }.epsilon(1e-10));
        CHECK(bases.zz == Approx{ D3 * D3 }.epsilon(1e-10));
        CHECK(bases.xy == 0);
        CHECK(bases.yz == 0);
        CHECK(bases.zx == 0);

        bases = mirror.covar_metric(curvi);
        CHECK(bases.xx == Approx{ D1 * D1 }.epsilon(1e-10));
        CHECK(bases.yy == Approx{ D2 * D2 }.epsilon(1e-10));
        CHECK(bases.zz == Approx{ D3 * D3 }.epsilon(1e-10));
        CHECK(bases.xy == 0);
        CHECK(bases.yz == 0);
        CHECK(bases.zx == 0);

        // contr metric
        bases = mirror.contr_metric(cart);
        CHECK(1 / bases.xx == Approx{ D1 * D1 }.epsilon(1e-10));
        CHECK(1 / bases.yy == Approx{ D2 * D2 }.epsilon(1e-10));
        CHECK(1 / bases.zz == Approx{ D3 * D3 }.epsilon(1e-10));
        CHECK(bases.xy == 0);
        CHECK(bases.yz == 0);
        CHECK(bases.zx == 0);

        bases = mirror.contr_metric(curvi);
        CHECK(1 / bases.xx == Approx{ D1 * D1 }.epsilon(1e-10));
        CHECK(1 / bases.yy == Approx{ D2 * D2 }.epsilon(1e-10));
        CHECK(1 / bases.zz == Approx{ D3 * D3 }.epsilon(1e-10));
        CHECK(bases.xy == 0);
        CHECK(bases.yz == 0);
        CHECK(bases.zx == 0);

        // covar basis
        basis = mirror.covar_basis<1>(cart);
        CHECK(basis.x == Approx{ D1 }.epsilon(1e-10));
        CHECK(basis.y == 0);
        CHECK(basis.z == 0);
        basis = mirror.covar_basis<2>(cart);
        CHECK(basis.x == 0);
        CHECK(basis.y == Approx{ D2 }.epsilon(1e-10));
        CHECK(basis.z == 0);
        basis = mirror.covar_basis<3>(cart);
        CHECK(basis.x == 0);
        CHECK(basis.y == 0);
        CHECK(basis.z == Approx{ D3 }.epsilon(1e-10));
        bases = mirror.covar_basis<0>(cart);
        CHECK(bases.xx == Approx{ D1 }.epsilon(1e-10));
        CHECK(bases.yy == Approx{ D2 }.epsilon(1e-10));
        CHECK(bases.zz == Approx{ D3 }.epsilon(1e-10));
        CHECK(bases.xy == 0);
        CHECK(bases.yz == 0);
        CHECK(bases.zx == 0);

        basis = mirror.covar_basis<1>(curvi);
        CHECK(basis.x == Approx{ D1 }.epsilon(1e-10));
        CHECK(basis.y == 0);
        CHECK(basis.z == 0);
        basis = mirror.covar_basis<2>(curvi);
        CHECK(basis.x == 0);
        CHECK(basis.y == Approx{ D2 }.epsilon(1e-10));
        CHECK(basis.z == 0);
        basis = mirror.covar_basis<3>(curvi);
        CHECK(basis.x == 0);
        CHECK(basis.y == 0);
        CHECK(basis.z == Approx{ D3 }.epsilon(1e-10));
        bases = mirror.covar_basis<0>(curvi);
        CHECK(bases.xx == Approx{ D1 }.epsilon(1e-10));
        CHECK(bases.yy == Approx{ D2 }.epsilon(1e-10));
        CHECK(bases.zz == Approx{ D3 }.epsilon(1e-10));
        CHECK(bases.xy == 0);
        CHECK(bases.yz == 0);
        CHECK(bases.zx == 0);

        // contr basis
        basis = mirror.contr_basis<1>(cart);
        CHECK(basis.x == Approx{ 1 / D1 }.epsilon(1e-10));
        CHECK(basis.y == 0);
        CHECK(basis.z == 0);
        basis = mirror.contr_basis<2>(cart);
        CHECK(basis.x == 0);
        CHECK(basis.y == Approx{ 1 / D2 }.epsilon(1e-10));
        CHECK(basis.z == 0);
        basis = mirror.contr_basis<3>(cart);
        CHECK(basis.x == 0);
        CHECK(basis.y == 0);
        CHECK(basis.z == Approx{ 1 / D3 }.epsilon(1e-10));
        bases = mirror.contr_basis<0>(cart);
        CHECK(bases.xx == Approx{ 1 / D1 }.epsilon(1e-10));
        CHECK(bases.yy == Approx{ 1 / D2 }.epsilon(1e-10));
        CHECK(bases.zz == Approx{ 1 / D3 }.epsilon(1e-10));
        CHECK(bases.xy == 0);
        CHECK(bases.yz == 0);
        CHECK(bases.zx == 0);

        basis = mirror.contr_basis<1>(curvi);
        CHECK(basis.x == Approx{ 1 / D1 }.epsilon(1e-10));
        CHECK(basis.y == 0);
        CHECK(basis.z == 0);
        basis = mirror.contr_basis<2>(curvi);
        CHECK(basis.x == 0);
        CHECK(basis.y == Approx{ 1 / D2 }.epsilon(1e-10));
        CHECK(basis.z == 0);
        basis = mirror.contr_basis<3>(curvi);
        CHECK(basis.x == 0);
        CHECK(basis.y == 0);
        CHECK(basis.z == Approx{ 1 / D3 }.epsilon(1e-10));
        bases = mirror.contr_basis<0>(curvi);
        CHECK(bases.xx == Approx{ 1 / D1 }.epsilon(1e-10));
        CHECK(bases.yy == Approx{ 1 / D2 }.epsilon(1e-10));
        CHECK(bases.zz == Approx{ 1 / D3 }.epsilon(1e-10));
        CHECK(bases.xy == 0);
        CHECK(bases.yz == 0);
        CHECK(bases.zx == 0);
    }

    { // inhomogeneous
        constexpr Real       xi = 0.512;
        constexpr Real       D1 = 2;
        constexpr Real       D2 = 0.43;
        constexpr Real       D3 = 1.54;
        MirrorGeometry const mirror{ xi, { D1, D2, D3 } };

        constexpr CartCoord cart{ 7.5121 };
        auto const          curvi = mirror.cotrans(cart);
        Vector              basis;
        Tensor              bases;

        // covar metric
        bases = mirror.covar_metric(cart);
        CHECK(bases.xx == Approx{ 997.702878094314 }.epsilon(1e-10));
        CHECK(bases.yy == Approx{ 0.011707557361683246 }.epsilon(1e-10));
        CHECK(bases.zz == Approx{ 0.15016572763097885 }.epsilon(1e-10));
        CHECK(bases.xy == 0);
        CHECK(bases.yz == 0);
        CHECK(bases.zx == 0);

        bases = mirror.covar_metric(curvi);
        CHECK(bases.xx == Approx{ 997.702878094314 }.epsilon(1e-10));
        CHECK(bases.yy == Approx{ 0.011707557361683246 }.epsilon(1e-10));
        CHECK(bases.zz == Approx{ 0.15016572763097885 }.epsilon(1e-10));
        CHECK(bases.xy == 0);
        CHECK(bases.yz == 0);
        CHECK(bases.zx == 0);

        // contr metric
        bases = mirror.contr_metric(cart);
        CHECK(bases.xx == Approx{ 0.0010023024108240257 }.epsilon(1e-10));
        CHECK(bases.yy == Approx{ 85.41491355599265 }.epsilon(1e-10));
        CHECK(bases.zz == Approx{ 6.659309123167078 }.epsilon(1e-10));
        CHECK(bases.xy == 0);
        CHECK(bases.yz == 0);
        CHECK(bases.zx == 0);

        bases = mirror.contr_metric(curvi);
        CHECK(bases.xx == Approx{ 0.0010023024108240257 }.epsilon(1e-10));
        CHECK(bases.yy == Approx{ 85.41491355599265 }.epsilon(1e-10));
        CHECK(bases.zz == Approx{ 6.659309123167078 }.epsilon(1e-10));
        CHECK(bases.xy == 0);
        CHECK(bases.yz == 0);
        CHECK(bases.zx == 0);

        // covar basis
        basis = mirror.covar_basis<1>(cart);
        CHECK(basis.x == Approx{ 31.586435033006083 }.epsilon(1e-10));
        CHECK(basis.y == 0);
        CHECK(basis.z == 0);
        basis = mirror.covar_basis<2>(cart);
        CHECK(basis.x == 0);
        CHECK(basis.y == Approx{ 0.10820146654127774 }.epsilon(1e-10));
        CHECK(basis.z == 0);
        basis = mirror.covar_basis<3>(cart);
        CHECK(basis.x == 0);
        CHECK(basis.y == 0);
        CHECK(basis.z == Approx{ 0.38751222900829707 }.epsilon(1e-10));
        bases = mirror.covar_basis<0>(cart);
        CHECK(bases.xx == Approx{ 31.586435033006083 }.epsilon(1e-10));
        CHECK(bases.yy == Approx{ 0.10820146654127774 }.epsilon(1e-10));
        CHECK(bases.zz == Approx{ 0.38751222900829707 }.epsilon(1e-10));
        CHECK(bases.xy == 0);
        CHECK(bases.yz == 0);
        CHECK(bases.zx == 0);

        basis = mirror.covar_basis<1>(curvi);
        CHECK(basis.x == Approx{ 31.586435033006083 }.epsilon(1e-10));
        CHECK(basis.y == 0);
        CHECK(basis.z == 0);
        basis = mirror.covar_basis<2>(curvi);
        CHECK(basis.x == 0);
        CHECK(basis.y == Approx{ 0.10820146654127774 }.epsilon(1e-10));
        CHECK(basis.z == 0);
        basis = mirror.covar_basis<3>(curvi);
        CHECK(basis.x == 0);
        CHECK(basis.y == 0);
        CHECK(basis.z == Approx{ 0.38751222900829707 }.epsilon(1e-10));
        bases = mirror.covar_basis<0>(curvi);
        CHECK(bases.xx == Approx{ 31.586435033006083 }.epsilon(1e-10));
        CHECK(bases.yy == Approx{ 0.10820146654127774 }.epsilon(1e-10));
        CHECK(bases.zz == Approx{ 0.38751222900829707 }.epsilon(1e-10));
        CHECK(bases.xy == 0);
        CHECK(bases.yz == 0);
        CHECK(bases.zx == 0);

        // contr basis
        basis = mirror.contr_basis<1>(cart);
        CHECK(basis.x == Approx{ 0.03165915998291846 }.epsilon(1e-10));
        CHECK(basis.y == 0);
        CHECK(basis.z == 0);
        basis = mirror.contr_basis<2>(cart);
        CHECK(basis.x == 0);
        CHECK(basis.y == Approx{ 9.24201891125487 }.epsilon(1e-10));
        CHECK(basis.z == 0);
        basis = mirror.contr_basis<3>(cart);
        CHECK(basis.x == 0);
        CHECK(basis.y == 0);
        CHECK(basis.z == Approx{ 2.5805637219737627 }.epsilon(1e-10));
        bases = mirror.contr_basis<0>(cart);
        CHECK(bases.xx == Approx{ 0.03165915998291846 }.epsilon(1e-10));
        CHECK(bases.yy == Approx{ 9.24201891125487 }.epsilon(1e-10));
        CHECK(bases.zz == Approx{ 2.5805637219737627 }.epsilon(1e-10));
        CHECK(bases.xy == 0);
        CHECK(bases.yz == 0);
        CHECK(bases.zx == 0);

        basis = mirror.contr_basis<1>(curvi);
        CHECK(basis.x == Approx{ 0.03165915998291846 }.epsilon(1e-10));
        CHECK(basis.y == 0);
        CHECK(basis.z == 0);
        basis = mirror.contr_basis<2>(curvi);
        CHECK(basis.x == 0);
        CHECK(basis.y == Approx{ 9.24201891125487 }.epsilon(1e-10));
        CHECK(basis.z == 0);
        basis = mirror.contr_basis<3>(curvi);
        CHECK(basis.x == 0);
        CHECK(basis.y == 0);
        CHECK(basis.z == Approx{ 2.5805637219737627 }.epsilon(1e-10));
        bases = mirror.contr_basis<0>(curvi);
        CHECK(bases.xx == Approx{ 0.03165915998291846 }.epsilon(1e-10));
        CHECK(bases.yy == Approx{ 9.24201891125487 }.epsilon(1e-10));
        CHECK(bases.zz == Approx{ 2.5805637219737627 }.epsilon(1e-10));
        CHECK(bases.xy == 0);
        CHECK(bases.yz == 0);
        CHECK(bases.zx == 0);
    }
}

TEST_CASE("Test LibPIC::MirrorGeometry::Transform", "[LibPIC::MirrorGeometry::Transform]")
{
    { // homogeneous
        constexpr Real       xi = 0;
        constexpr Real       D1 = 0.1;
        constexpr Real       D2 = 0.43;
        constexpr Real       D3 = 1.54;
        MirrorGeometry const mirror{ xi, { D1, D2, D3 } };

        constexpr CartCoord cart{ 14.5 };
        auto const          curvi = mirror.cotrans(cart);
        constexpr Vector    vec{ 1.3, .506, -.598 };
        Vector              tmp;

        tmp = mirror.contr_to_covar(vec, cart);
        CHECK(tmp.x == Approx{ vec.x * D1 * D1 }.epsilon(1e-10));
        CHECK(tmp.y == Approx{ vec.y * D2 * D2 }.epsilon(1e-10));
        CHECK(tmp.z == Approx{ vec.z * D3 * D3 }.epsilon(1e-10));
        tmp = mirror.contr_to_covar(vec, curvi);
        CHECK(tmp.x == Approx{ vec.x * D1 * D1 }.epsilon(1e-10));
        CHECK(tmp.y == Approx{ vec.y * D2 * D2 }.epsilon(1e-10));
        CHECK(tmp.z == Approx{ vec.z * D3 * D3 }.epsilon(1e-10));

        tmp = mirror.covar_to_contr(vec, cart);
        CHECK(tmp.x == Approx{ vec.x / D1 / D1 }.epsilon(1e-10));
        CHECK(tmp.y == Approx{ vec.y / D2 / D2 }.epsilon(1e-10));
        CHECK(tmp.z == Approx{ vec.z / D3 / D3 }.epsilon(1e-10));
        tmp = mirror.covar_to_contr(vec, curvi);
        CHECK(tmp.x == Approx{ vec.x / D1 / D1 }.epsilon(1e-10));
        CHECK(tmp.y == Approx{ vec.y / D2 / D2 }.epsilon(1e-10));
        CHECK(tmp.z == Approx{ vec.z / D3 / D3 }.epsilon(1e-10));

        tmp = mirror.cart_to_contr(vec, cart);
        CHECK(tmp.x == Approx{ vec.x / D1 }.epsilon(1e-10));
        CHECK(tmp.y == Approx{ vec.y / D2 }.epsilon(1e-10));
        CHECK(tmp.z == Approx{ vec.z / D3 }.epsilon(1e-10));
        tmp = mirror.cart_to_contr(vec, curvi);
        CHECK(tmp.x == Approx{ vec.x / D1 }.epsilon(1e-10));
        CHECK(tmp.y == Approx{ vec.y / D2 }.epsilon(1e-10));
        CHECK(tmp.z == Approx{ vec.z / D3 }.epsilon(1e-10));

        tmp = mirror.contr_to_cart(vec, cart);
        CHECK(tmp.x == Approx{ vec.x * D1 }.epsilon(1e-10));
        CHECK(tmp.y == Approx{ vec.y * D2 }.epsilon(1e-10));
        CHECK(tmp.z == Approx{ vec.z * D3 }.epsilon(1e-10));
        tmp = mirror.contr_to_cart(vec, curvi);
        CHECK(tmp.x == Approx{ vec.x * D1 }.epsilon(1e-10));
        CHECK(tmp.y == Approx{ vec.y * D2 }.epsilon(1e-10));
        CHECK(tmp.z == Approx{ vec.z * D3 }.epsilon(1e-10));

        tmp = mirror.cart_to_covar(vec, cart);
        CHECK(tmp.x == Approx{ vec.x * D1 }.epsilon(1e-10));
        CHECK(tmp.y == Approx{ vec.y * D2 }.epsilon(1e-10));
        CHECK(tmp.z == Approx{ vec.z * D3 }.epsilon(1e-10));
        tmp = mirror.cart_to_covar(vec, curvi);
        CHECK(tmp.x == Approx{ vec.x * D1 }.epsilon(1e-10));
        CHECK(tmp.y == Approx{ vec.y * D2 }.epsilon(1e-10));
        CHECK(tmp.z == Approx{ vec.z * D3 }.epsilon(1e-10));

        tmp = mirror.covar_to_cart(vec, cart);
        CHECK(tmp.x == Approx{ vec.x / D1 }.epsilon(1e-10));
        CHECK(tmp.y == Approx{ vec.y / D2 }.epsilon(1e-10));
        CHECK(tmp.z == Approx{ vec.z / D3 }.epsilon(1e-10));
        tmp = mirror.covar_to_cart(vec, curvi);
        CHECK(tmp.x == Approx{ vec.x / D1 }.epsilon(1e-10));
        CHECK(tmp.y == Approx{ vec.y / D2 }.epsilon(1e-10));
        CHECK(tmp.z == Approx{ vec.z / D3 }.epsilon(1e-10));
    }

    { // inhomogeneous
        constexpr Real       xi = 0.512;
        constexpr Real       D1 = 2;
        constexpr Real       D2 = 0.43;
        constexpr Real       D3 = 1.54;
        MirrorGeometry const mirror{ xi, { D1, D2, D3 } };

        constexpr CartCoord cart{ 7.5121 };
        auto const          curvi = mirror.cotrans(cart);
        constexpr Vector    vec{ 1.3, .506, -.598 };
        Vector              tmp;

        tmp = mirror.contr_to_covar(vec, cart);
        CHECK(tmp.x == Approx{ 1297.0137415226081 }.epsilon(1e-10));
        CHECK(tmp.y == Approx{ 0.005924024025011723 }.epsilon(1e-10));
        CHECK(tmp.z == Approx{ -0.08979910512332535 }.epsilon(1e-10));
        tmp = mirror.contr_to_covar(vec, curvi);
        CHECK(tmp.x == Approx{ 1297.0137415226081 }.epsilon(1e-10));
        CHECK(tmp.y == Approx{ 0.005924024025011723 }.epsilon(1e-10));
        CHECK(tmp.z == Approx{ -0.08979910512332535 }.epsilon(1e-10));

        tmp = mirror.covar_to_contr(vec, cart);
        CHECK(tmp.x == Approx{ 0.0013029931340712334 }.epsilon(1e-10));
        CHECK(tmp.y == Approx{ 43.21994625933228 }.epsilon(1e-10));
        CHECK(tmp.z == Approx{ -3.9822668556539123 }.epsilon(1e-10));
        tmp = mirror.covar_to_contr(vec, curvi);
        CHECK(tmp.x == Approx{ 0.0013029931340712334 }.epsilon(1e-10));
        CHECK(tmp.y == Approx{ 43.21994625933228 }.epsilon(1e-10));
        CHECK(tmp.z == Approx{ -3.9822668556539123 }.epsilon(1e-10));

        tmp = mirror.cart_to_contr(vec, cart);
        CHECK(tmp.x == Approx{ 0.041156907977794 }.epsilon(1e-10));
        CHECK(tmp.y == Approx{ 4.676461569094965 }.epsilon(1e-10));
        CHECK(tmp.z == Approx{ -1.54317710574031 }.epsilon(1e-10));
        tmp = mirror.cart_to_contr(vec, curvi);
        CHECK(tmp.x == Approx{ 0.041156907977794 }.epsilon(1e-10));
        CHECK(tmp.y == Approx{ 4.676461569094965 }.epsilon(1e-10));
        CHECK(tmp.z == Approx{ -1.54317710574031 }.epsilon(1e-10));

        tmp = mirror.contr_to_cart(vec, cart);
        CHECK(tmp.x == Approx{ 41.06236554290791 }.epsilon(1e-10));
        CHECK(tmp.y == Approx{ 0.05474994206988654 }.epsilon(1e-10));
        CHECK(tmp.z == Approx{ -0.23173231294696164 }.epsilon(1e-10));
        tmp = mirror.contr_to_cart(vec, curvi);
        CHECK(tmp.x == Approx{ 41.06236554290791 }.epsilon(1e-10));
        CHECK(tmp.y == Approx{ 0.05474994206988654 }.epsilon(1e-10));
        CHECK(tmp.z == Approx{ -0.23173231294696164 }.epsilon(1e-10));

        tmp = mirror.cart_to_covar(vec, cart);
        CHECK(tmp.x == Approx{ 41.06236554290791 }.epsilon(1e-10));
        CHECK(tmp.y == Approx{ 0.05474994206988654 }.epsilon(1e-10));
        CHECK(tmp.z == Approx{ -0.23173231294696164 }.epsilon(1e-10));
        tmp = mirror.cart_to_covar(vec, curvi);
        CHECK(tmp.x == Approx{ 41.06236554290791 }.epsilon(1e-10));
        CHECK(tmp.y == Approx{ 0.05474994206988654 }.epsilon(1e-10));
        CHECK(tmp.z == Approx{ -0.23173231294696164 }.epsilon(1e-10));

        tmp = mirror.covar_to_cart(vec, cart);
        CHECK(tmp.x == Approx{ 0.041156907977794 }.epsilon(1e-10));
        CHECK(tmp.y == Approx{ 4.676461569094965 }.epsilon(1e-10));
        CHECK(tmp.z == Approx{ -1.54317710574031 }.epsilon(1e-10));
        tmp = mirror.covar_to_cart(vec, curvi);
        CHECK(tmp.x == Approx{ 0.041156907977794 }.epsilon(1e-10));
        CHECK(tmp.y == Approx{ 4.676461569094965 }.epsilon(1e-10));
        CHECK(tmp.z == Approx{ -1.54317710574031 }.epsilon(1e-10));
    }
}
