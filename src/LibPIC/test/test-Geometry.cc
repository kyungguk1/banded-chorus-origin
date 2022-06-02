/*
 * Copyright (c) 2021-2022, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "catch2/catch.hpp"

#define LIBPIC_INLINE_VERSION 1
#include <PIC/Geometry.h>

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
