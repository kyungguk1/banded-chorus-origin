/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "catch2/catch.hpp"

#include <PIC/FourTensor.h>
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

TEST_CASE("Test libPIC::FourTensor", "[libPIC::FourTensor]")
{
    {
        constexpr FourTensor t1{};
        CHECK(*t1.tt == 0);
        CHECK((t1.ts.x == 0 && t1.ts.y == 0 && t1.ts.z == 0));
        CHECK((t1.ss.xx == 0 && t1.ss.yy == 0 && t1.ss.zz == 0));
        CHECK((t1.ss.xy == 0 && t1.ss.yz == 0 && t1.ss.zx == 0));

        constexpr FourTensor t2{ 1 };
        CHECK(*t2.tt == 1);
        CHECK((t2.ts.x == 1 && t2.ts.y == 1 && t2.ts.z == 1));
        CHECK((t2.ss.xx == 1 && t2.ss.yy == 1 && t2.ss.zz == 1));
        CHECK((t2.ss.xy == 1 && t2.ss.yz == 1 && t2.ss.zx == 1));

        constexpr FourTensor t3{ -1, { 10, 20, 30 }, { 1, 2, 3, 4, 5, 6 } };
        CHECK(*t3.tt == -1);
        CHECK((t3.ts.x == 10 && t3.ts.y == 20 && t3.ts.z == 30));
        CHECK((t3.ss.xx == 1 && t3.ss.yy == 2 && t3.ss.zz == 3));
        CHECK((t3.ss.xy == 4 && t3.ss.yz == 5 && t3.ss.zx == 6));

        constexpr bool tf2 = std::addressof(t1) == std::addressof(+t1);
        CHECK(tf2);

        constexpr FourTensor t4 = -t3;
        CHECK(*t4.tt == 1);
        CHECK((t4.ts.x == -10 && t4.ts.y == -20 && t4.ts.z == -30));
        CHECK((t4.ss.xx == -1 && t4.ss.yy == -2 && t4.ss.zz == -3));
        CHECK((t4.ss.xy == -4 && t4.ss.yz == -5 && t4.ss.zx == -6));
    }

    {
        constexpr auto is_equal = [](FourTensor lhs, FourTensor rhs) {
            return *lhs.tt == *rhs.tt && lhs.ts.x == rhs.ts.x && lhs.ts.y == rhs.ts.y
                && lhs.ts.z == rhs.ts.z && lhs.ss.xx == rhs.ss.xx && lhs.ss.yy == rhs.ss.yy
                && lhs.ss.zz == rhs.ss.zz && lhs.ss.xy == rhs.ss.xy && lhs.ss.yz == rhs.ss.yz
                && lhs.ss.zx == rhs.ss.zx;
        };

        constexpr FourTensor v1{ -1, { 1, 2, 3 }, { 1, 2, 3, -1, -2, -3 } };
        constexpr double     x{ 1 };
        CHECK(is_equal(v1 + x, { 0, { 2, 3, 4 }, { 2, 3, 4, -0, -1, -2 } }));
        CHECK(is_equal(v1 - x, { -2, { 0, 1, 2 }, { 0, 1, 2, -2, -3, -4 } }));
        CHECK(is_equal(v1 * x, v1));
        CHECK(is_equal(v1 / x, v1));
        CHECK(is_equal(x + v1, v1 + x));
        CHECK(is_equal(x - v1, -(v1 - x)));
        CHECK(is_equal(x * v1, v1));
        CHECK(is_equal(
            x / v1,
            { -x / 1, { x / 1, x / 2, x / 3 }, { x / 1, x / 2, x / 3, -x / 1, -x / 2, -x / 3 } }));

        constexpr FourTensor v2 = v1 * 10.;
        CHECK(is_equal(v1 + v2, { -11, { 11, 22, 33 }, { 11, 22, 33, -11, -22, -33 } }));
        CHECK(is_equal(v2 - v1, { -9, { 9, 18, 27 }, { 9, 18, 27, -9, -18, -27 } }));
        CHECK(is_equal(v1 * v2, { 10, { 10, 40, 90 }, { 10, 40, 90, 10, 40, 90 } }));
        CHECK(is_equal(v2 / v1, { 10, { 10, 10, 10 }, { 10, 10, 10, 10, 10, 10 } }));
    }
}