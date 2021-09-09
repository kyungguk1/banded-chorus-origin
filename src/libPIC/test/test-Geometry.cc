/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "catch2/catch.hpp"

#include <PIC/Geometry.h>

TEST_CASE("Test libPIC::Geometry", "[libPIC::Geometry]")
{
    REQUIRE_THROWS_AS(Geometry({ 0, 0, 1 }), std::exception);
    REQUIRE_NOTHROW(Geometry({ 1, 0, 0 }));
    REQUIRE(Geometry::e3.x == 0);
    REQUIRE(Geometry::e3.y == 0);
    REQUIRE(Geometry::e3.z == 1);

    constexpr Vector B0{ 2, 5, 0 };
    auto const       mag   = std::sqrt(dot(B0, B0));
    auto const       e1    = B0 / mag;
    auto const       theta = std::atan2(e1.y, e1.x);
    Geometry const   geo{ mag, theta };

    REQUIRE(dot(geo.B0, geo.B0) == Approx{ dot(B0, B0) }.epsilon(1e-15));
    REQUIRE(dot(geo.e1, geo.e1) == Approx{ dot(e1, e1) }.epsilon(1e-15));

    auto const e3 = cross(geo.e1, geo.e2);
    REQUIRE(dot(geo.e3, geo.e3) == Approx{ dot(e3, e3) }.epsilon(1e-15));

    auto const cart_e1 = geo.fac2cart({ 1, 0, 0 });
    CHECK(dot(cart_e1, cart_e1) == Approx{ dot(e1, e1) }.epsilon(1e-15));
    auto const delta_fac_e1 = geo.cart2fac(e1) - Vector{ 1, 0, 0 };
    CHECK(dot(delta_fac_e1, delta_fac_e1) == Approx{ 0 }.margin(1e-15));

    Vector const vfac{ 0.30, 0.5604, 0.493 };
    Tensor const tfac{ 0.30, 0.5604, 0.493, 0, 0, 0 };

    auto const vcart = geo.fac2cart(vfac);
    auto const tcart = geo.fac2cart(tfac);
    CHECK(dot(vcart, vcart) == Approx{ dot(vfac, vfac) }.epsilon(1e-15));
    CHECK(trace(tfac) == Approx{ trace(tcart) }.epsilon(1e-15));

    auto const d_vfac = geo.cart2fac(vcart) - vfac;
    auto const d_tfac = geo.cart2fac(tcart) - tfac;
    CHECK(dot(d_vfac, d_vfac) == Approx{ 0 }.margin(1e-15));
    CHECK(dot(d_tfac.lo(), d_tfac.lo()) == Approx{ 0 }.margin(1e-15));
    CHECK(dot(d_tfac.hi(), d_tfac.hi()) == Approx{ 0 }.margin(1e-15));

    FourVector const Vfac{ 1, vfac };
    FourTensor const Tfac{ 1, vfac, tfac };

    auto const Vcart = geo.fac2cart(Vfac);
    CHECK(*Vcart.t == *Vfac.t);
    CHECK(Vcart.s.x == vcart.x);
    CHECK(Vcart.s.y == vcart.y);
    CHECK(Vcart.s.z == vcart.z);
    auto const Tcart = geo.fac2cart(Tfac);
    CHECK(*Tcart.tt == *Tfac.tt);
    CHECK(Tcart.ts.x == vcart.x);
    CHECK(Tcart.ts.y == vcart.y);
    CHECK(Tcart.ts.z == vcart.z);
    CHECK(Tcart.ss.xx == tcart.xx);
    CHECK(Tcart.ss.yy == tcart.yy);
    CHECK(Tcart.ss.zz == tcart.zz);
    CHECK(Tcart.ss.xy == tcart.xy);
    CHECK(Tcart.ss.yz == tcart.yz);
    CHECK(Tcart.ss.zx == tcart.zx);

    auto const Vfac2 = geo.cart2fac(Vcart);
    CHECK(*Vfac2.t == *Vfac.t);
    CHECK(Vfac2.s.x == Approx{ Vfac.s.x }.epsilon(1e-15));
    CHECK(Vfac2.s.y == Approx{ Vfac.s.y }.epsilon(1e-15));
    CHECK(Vfac2.s.z == Approx{ Vfac.s.z }.epsilon(1e-15));
    auto const Tfac2 = geo.cart2fac(Tcart);
    CHECK(*Tfac2.tt == *Tfac.tt);
    CHECK(Tfac2.ts.x == Approx{ Tfac.ts.x }.epsilon(1e-15));
    CHECK(Tfac2.ts.y == Approx{ Tfac.ts.y }.epsilon(1e-15));
    CHECK(Tfac2.ts.z == Approx{ Tfac.ts.z }.epsilon(1e-15));
    CHECK(Tfac2.ss.xx == Approx{ Tfac.ss.xx }.epsilon(1e-15));
    CHECK(Tfac2.ss.yy == Approx{ Tfac.ss.yy }.epsilon(1e-15));
    CHECK(Tfac2.ss.zz == Approx{ Tfac.ss.zz }.epsilon(1e-15));
    CHECK(Tfac2.ss.xy == Approx{ Tfac.ss.xy }.epsilon(1e-15));
    CHECK(Tfac2.ss.yz == Approx{ Tfac.ss.yz }.epsilon(1e-15));
    CHECK(Tfac2.ss.zx == Approx{ Tfac.ss.zx }.epsilon(1e-15));
}
