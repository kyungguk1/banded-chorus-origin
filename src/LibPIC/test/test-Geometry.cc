/*
 * Copyright (c) 2021-2022, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "catch2/catch.hpp"

#define LIBPIC_INLINE_VERSION 1
#include <PIC/Geometry.h>

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

TEST_CASE("Test LibPIC::Geometry", "[LibPIC::Geometry]")
{
    constexpr Real xi = 0.512;
    constexpr Real D1 = 2;
    constexpr Real B0 = 2.40939;
    Geometry const geo{ xi, D1, B0 };

    CHECK(geo.B0() == B0);

    Vector     v;
    Tensor     t;
    FourVector fv;
    FourTensor ft;
    CHECK(&v == &geo.mfa_to_cart(v, CartCoord{}));
    CHECK(&v == &geo.mfa_to_cart(v, CurviCoord{}));
    CHECK(&t == &geo.mfa_to_cart(t, CartCoord{}));
    CHECK(&t == &geo.mfa_to_cart(t, CurviCoord{}));
    CHECK(&fv == &geo.mfa_to_cart(fv, CartCoord{}));
    CHECK(&fv == &geo.mfa_to_cart(fv, CurviCoord{}));
    CHECK(&ft == &geo.mfa_to_cart(ft, CartCoord{}));
    CHECK(&ft == &geo.mfa_to_cart(ft, CurviCoord{}));

    CHECK(&v == &geo.cart_to_mfa(v, CartCoord{}));
    CHECK(&v == &geo.cart_to_mfa(v, CurviCoord{}));
    CHECK(&t == &geo.cart_to_mfa(t, CartCoord{}));
    CHECK(&t == &geo.cart_to_mfa(t, CurviCoord{}));
    CHECK(&fv == &geo.cart_to_mfa(fv, CartCoord{}));
    CHECK(&fv == &geo.cart_to_mfa(fv, CurviCoord{}));
    CHECK(&ft == &geo.cart_to_mfa(ft, CartCoord{}));
    CHECK(&ft == &geo.cart_to_mfa(ft, CurviCoord{}));

    CHECK(geo.e1(CurviCoord{ 0 }).x == 1);
    CHECK(geo.e1(CurviCoord{ 0 }).y == 0);
    CHECK(geo.e1(CurviCoord{ 0 }).z == 0);
    CHECK(geo.e2(CurviCoord{ 0 }).x == 0);
    CHECK(geo.e2(CurviCoord{ 0 }).y == 1);
    CHECK(geo.e2(CurviCoord{ 0 }).z == 0);
    CHECK(geo.e3(CurviCoord{ 0 }).x == 0);
    CHECK(geo.e3(CurviCoord{ 0 }).y == 0);
    CHECK(geo.e3(CurviCoord{ 0 }).z == 1);

    constexpr CartCoord cart{ 4.5121 };
    auto const          curvi = geo.cotrans(cart);
    Vector              B1;
    Vector              B2;

    B1 = geo.Bcontr(cart);
    CHECK(B1.x / B0 == Approx{ geo.Bcontr_div_B0(curvi).x }.epsilon(1e-10));
    CHECK(B1.y == 0);
    CHECK(B1.z == 0);
    B1 = geo.Bcontr(curvi);
    CHECK(B1.x / B0 == Approx{ geo.Bcontr_div_B0(cart).x }.epsilon(1e-10));
    CHECK(B1.y == 0);
    CHECK(B1.z == 0);

    B1 = geo.Bcovar(cart);
    CHECK(B1.x / B0 == Approx{ geo.Bcovar_div_B0(curvi).x }.epsilon(1e-10));
    CHECK(B1.y == 0);
    CHECK(B1.z == 0);
    B1 = geo.Bcovar(curvi);
    CHECK(B1.x / B0 == Approx{ geo.Bcovar_div_B0(cart).x }.epsilon(1e-10));
    CHECK(B1.y == 0);
    CHECK(B1.z == 0);

    B1 = geo.Bcart(cart);
    CHECK(B1.x / B0 == Approx{ geo.Bcart_div_B0(curvi).x }.epsilon(1e-10));
    CHECK(B1.y == 0);
    CHECK(B1.z == 0);
    B1 = geo.Bcart(curvi);
    CHECK(B1.x / B0 == Approx{ geo.Bcart_div_B0(cart).x }.epsilon(1e-10));
    CHECK(B1.y == 0);
    CHECK(B1.z == 0);

    CHECK(std::pow(geo.Bmag(cart), 2) == Approx{ dot(B1, B1) }.epsilon(1e-10));
    CHECK(std::pow(geo.Bmag(curvi), 2) == Approx{ dot(B1, B1) }.epsilon(1e-10));

    B1 = geo.Bcart(cart, 3, 2);
    B2 = geo.Bcart_div_B0(curvi, 3, 2);
    CHECK(B1.x / B0 == Approx{ B2.x }.epsilon(1e-10));
    CHECK(B1.y / B0 == Approx{ B2.y }.epsilon(1e-10));
    CHECK(B1.z / B0 == Approx{ B2.z }.epsilon(1e-10));
    B1 = geo.Bcart(curvi, 3, 2);
    B2 = geo.Bcart_div_B0(cart, 3, 2);
    CHECK(B1.x / B0 == Approx{ B2.x }.epsilon(1e-10));
    CHECK(B1.y / B0 == Approx{ B2.y }.epsilon(1e-10));
    CHECK(B1.z / B0 == Approx{ B2.z }.epsilon(1e-10));
}

#if 0
#include <PIC/BorisPush.h>
#include <PIC/println.h>
#include <fstream>
#include <vector>

TEST_CASE("Test LibPIC::Geometry::MirrorMotion", "[LibPIC::Geometry::MirrorMotion]")
{
    constexpr Real  O0 = M_PI * 2, ob0 = O0 * 0.05;
    constexpr Real  x0 = 0, v0 = M_PI * 2, pa0 = 70 * M_PI / 180, ph0 = 0 * M_PI / 180;
    constexpr Real  xi = ob0 / v0, D1 = 1;
    constexpr Real  Dt = M_PI * 2 / O0 / 100;
    constexpr Real  nt = (M_PI * 2 / O0 / Dt) * (2 * O0 / ob0);
    Geometry const  geo{ xi, D1, O0 };
    BorisPush const boris{ Dt, 1, O0, O0 };

    Vector vv = { v0 * std::cos(pa0), v0 * std::sin(pa0) * std::cos(ph0), v0 * std::sin(pa0) * std::sin(ph0) };
    Real   x  = x0;

    std::vector<decltype(vv)> vel{ vv };
    std::vector<decltype(x)>  pos{ x };
    for (int i = 0; i < nt; ++i) {
        x += .5 * vv.x * Dt;
        {
            auto xx = Vector{ x, 0, 0 };
            xx      = xx - cross(vv, { 1, 0, 0 }) / (O0 * (1 + xi * xi * x * x));
            boris.non_relativistic(vv, geo.Bcart(CartCoord{ x }, xx.y, xx.z), { 0, 0, 0 });
        }
        x += .5 * vv.x * Dt;

        vel.push_back(vv);
        pos.push_back(x);
    }

    const auto printer = [&] {
        std::ofstream os{ "/Users/kyungguk/Downloads/mirror_motion.m" };
        println(os, '{');
        for (unsigned i = 0; i < vel.size() - 1; ++i) {
            println(os, "    {", vel[i], ", ", pos[i], "}, ");
        }
        println(os, "    {", vel.back(), ", ", pos.back(), '}');
        println(os, '}');
        os.close();
    };
    printer();
}
#endif
