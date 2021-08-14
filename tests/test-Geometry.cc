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

    REQUIRE(std::abs(dot(geo.B0 - B0, geo.B0 - B0)) < mag * mag * 1e-15);
    REQUIRE(std::abs(dot(geo.e1 - e1, geo.e1 - e1)) < 1e-15);

    auto const e3 = cross(geo.e1, geo.e2);
    REQUIRE(std::abs(dot(geo.e3 - e3, geo.e3 - e3)) < 1e-15);

    auto const cart_e1 = geo.fac2cart({ 1, 0, 0 });
    CHECK(std::abs(dot(cart_e1 - e1, cart_e1 - e1)) < 1e-15);
    auto const delta_fac_e1 = geo.cart2fac(e1) - Vector{ 1, 0, 0 };
    CHECK(std::abs(dot(delta_fac_e1, delta_fac_e1)) < 1e-15);

    Vector const vfac{ 0.30, 0.5604, 0.493 };
    Tensor const tfac{ 0.30, 0.5604, 0.493, 0, 0, 0 };

    auto const vcart = geo.fac2cart(vfac);
    auto const tcart = geo.fac2cart(tfac);
    CHECK(std::abs(dot(vcart, vcart) - dot(vfac, vfac)) < 1e-15);
    CHECK(std::abs(trace(tfac) - trace(tcart)) < 1e-15);

    auto const d_vfac = geo.cart2fac(vcart) - vfac;
    auto const d_tfac = geo.cart2fac(tcart) - tfac;
    CHECK(std::abs(dot(d_vfac, d_vfac)) < 1e-15);
    CHECK(std::abs(dot(d_tfac.lo(), d_tfac.lo())) < 1e-15);
    CHECK(std::abs(dot(d_tfac.hi(), d_tfac.hi())) < 1e-15);
}
