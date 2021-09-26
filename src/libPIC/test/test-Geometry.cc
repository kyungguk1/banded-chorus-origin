/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "catch2/catch.hpp"

#include <PIC/Geometry.h>

TEST_CASE("Test libPIC::Geometry", "[libPIC::Geometry]")
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
    CHECK(&v == &geo.fac_to_cart(v, CartCoord{}));
    CHECK(&v == &geo.fac_to_cart(v, CurviCoord{}));
    CHECK(&t == &geo.fac_to_cart(t, CartCoord{}));
    CHECK(&t == &geo.fac_to_cart(t, CurviCoord{}));
    CHECK(&fv == &geo.fac_to_cart(fv, CartCoord{}));
    CHECK(&fv == &geo.fac_to_cart(fv, CurviCoord{}));
    CHECK(&ft == &geo.fac_to_cart(ft, CartCoord{}));
    CHECK(&ft == &geo.fac_to_cart(ft, CurviCoord{}));

    CHECK(&v == &geo.cart_to_fac(v, CartCoord{}));
    CHECK(&v == &geo.cart_to_fac(v, CurviCoord{}));
    CHECK(&t == &geo.cart_to_fac(t, CartCoord{}));
    CHECK(&t == &geo.cart_to_fac(t, CurviCoord{}));
    CHECK(&fv == &geo.cart_to_fac(fv, CartCoord{}));
    CHECK(&fv == &geo.cart_to_fac(fv, CurviCoord{}));
    CHECK(&ft == &geo.cart_to_fac(ft, CartCoord{}));
    CHECK(&ft == &geo.cart_to_fac(ft, CurviCoord{}));

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
