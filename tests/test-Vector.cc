/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "catch2/catch.hpp"

#include <Utility/Vector.h>

TEST_CASE("Test libPIC::Vector", "[libPIC::Vector]")
{
    {
        constexpr Vector v1{};
        constexpr bool   tf = v1.fold(true, [](bool lhs, auto rhs) {
            return lhs && rhs == 0;
          });
        CHECK(tf);

        constexpr Vector v2{ 1 };
        CHECK(v2.fold(true, [](bool lhs, auto rhs) {
            return lhs && rhs == 1;
        }));

        constexpr Vector v3{ 1, 2, 3 };
        CHECK((v3.x == 1 && v3.y == 2 && v3.z == 3));

        constexpr bool tf2 = std::addressof(v1) == std::addressof(+v1);
        CHECK(tf2);

        constexpr Vector v4 = -v3;
        CHECK((v4.x == -1 && v4.y == -2 && v4.z == -3));

        constexpr auto dot1 = dot(v1, v3);
        CHECK(dot1 == 0);
        constexpr auto dot2 = dot(v2, v3);
        CHECK(dot2 == 6);

        constexpr auto cross1 = cross(v1, v3);
        CHECK(cross1.fold(true, [](bool lhs, auto rhs) {
            return lhs && rhs == 0;
        }));
        constexpr auto cross2 = cross(Vector{ 4, 3, 8 }, v3);
        CHECK((cross2.x == -7 && cross2.y == -4 && cross2.z == 5));
    }

    {
        constexpr auto is_equal = [](Vector lhs, Vector rhs) {
            return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
        };

        constexpr Vector v1{ 1, 2, 3 };
        constexpr double x{ 1 };
        CHECK(is_equal(v1 + x, { 2, 3, 4 }));
        CHECK(is_equal(v1 - x, { 0, 1, 2 }));
        CHECK(is_equal(v1 * x, v1));
        CHECK(is_equal(v1 / x, v1));
        CHECK(is_equal(x + v1, v1 + x));
        CHECK(is_equal(x - v1, -(v1 - x)));
        CHECK(is_equal(x * v1, v1));
        CHECK(is_equal(x / v1, { x / 1, x / 2, x / 3 }));

        constexpr Vector v2 = v1 * 10.;
        CHECK(is_equal(v1 + v2, { 11, 22, 33 }));
        CHECK(is_equal(v2 - v1, { 9, 18, 27 }));
        CHECK(is_equal(v1 * v2, { 10, 40, 90 }));
        CHECK(is_equal(v2 / v1, { 10, 10, 10 }));
    }
}
