/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "catch2/catch.hpp"

#include <PIC/Tensor.h>

TEST_CASE("Test libPIC::Tensor", "[libPIC::Tensor]")
{
    {
        constexpr Tensor t1{};
        constexpr bool   tf = t1.fold(true, [](bool lhs, auto rhs) {
            return lhs && rhs == 0;
          });
        CHECK(tf);

        constexpr Tensor t2{ 1 };
        CHECK(t2.fold(true, [](bool lhs, auto rhs) {
            return lhs && rhs == 1;
        }));

        constexpr Tensor t3{ 1, 2, 3, 4, 5, 6 };
        CHECK((t3.xx == 1 && t3.yy == 2 && t3.zz == 3));
        CHECK((t3.xy == 4 && t3.yz == 5 && t3.zx == 6));

        constexpr bool tf2 = std::addressof(t1) == std::addressof(+t1);
        CHECK(tf2);

        constexpr Tensor t4 = -t3;
        CHECK((t4.xx == -1 && t4.yy == -2 && t4.zz == -3));
        CHECK((t4.xy == -4 && t4.yz == -5 && t4.zx == -6));

        constexpr auto is_equal = [](Vector lhs, Vector rhs) {
            return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
        };
        constexpr Vector v1{ 2, 4, 5 };
        constexpr auto   dot1 = dot(v1, t3);
        CHECK(is_equal(dot1, { 48, 41, 47 }));
        constexpr auto dot2 = dot(t3, v1);
        CHECK(is_equal(dot2, { 48, 41, 47 }));

        CHECK(is_equal(t3.lo(), { 1, 2, 3 }));
        CHECK(is_equal(t3.hi(), { 4, 5, 6 }));

        constexpr auto tr = trace(t3);
        CHECK(std::abs(tr - t3.lo().fold(double{}, std::plus{})) < 1e-15);
    }

    {
        constexpr auto is_equal = [](Tensor lhs, Tensor rhs) {
            return lhs.xx == rhs.xx && lhs.yy == rhs.yy && lhs.zz == rhs.zz && lhs.xy == rhs.xy
                && lhs.yz == rhs.yz && lhs.zx == rhs.zx;
        };

        constexpr Tensor v1{ 1, 2, 3, -1, -2, -3 };
        constexpr double x{ 1 };
        CHECK(is_equal(v1 + x, { 2, 3, 4, -0, -1, -2 }));
        CHECK(is_equal(v1 - x, { 0, 1, 2, -2, -3, -4 }));
        CHECK(is_equal(v1 * x, v1));
        CHECK(is_equal(v1 / x, v1));
        CHECK(is_equal(x + v1, v1 + x));
        CHECK(is_equal(x - v1, -(v1 - x)));
        CHECK(is_equal(x * v1, v1));
        CHECK(is_equal(x / v1, { x / 1, x / 2, x / 3, -x / 1, -x / 2, -x / 3 }));

        constexpr Tensor v2 = v1 * 10.;
        CHECK(is_equal(v1 + v2, { 11, 22, 33, -11, -22, -33 }));
        CHECK(is_equal(v2 - v1, { 9, 18, 27, -9, -18, -27 }));
        CHECK(is_equal(v1 * v2, { 10, 40, 90, 10, 40, 90 }));
        CHECK(is_equal(v2 / v1, { 10, 10, 10, 10, 10, 10 }));
    }
}
