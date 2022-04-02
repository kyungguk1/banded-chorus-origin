/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "catch2/catch.hpp"

#include <PIC/CartCoord.h>
#include <PIC/CurviCoord.h>
#include <PIC/Scalar.h>

TEST_CASE("Test LibPIC::Scalar", "[LibPIC::Scalar]")
{
    {
        constexpr Scalar s1{};
        CHECK(0 == *s1);
        CHECK(*s1 == s1());
        CHECK(*s1 == static_cast<double>(s1));

        constexpr Scalar s2{ 1 };
        CHECK(1 == *s2);
        CHECK(*s2 == s2());
        CHECK(*s2 == static_cast<double>(s2));

        constexpr auto s3 = -s2;
        CHECK(-1 == *s3);
        CHECK(*s3 == s3());
        CHECK(*s3 == static_cast<double>(s3));
    }

    {
        constexpr Scalar s1{ 1 };
        constexpr long   s2{ 2 };
        CHECK(3 == *(s1 + s2));
        CHECK(-1 == *(s1 - s2));
        CHECK(2 == *(s1 * s2));
        CHECK(.5 == *(s1 / s2));
    }
    {
        constexpr long   s1{ 1 };
        constexpr Scalar s2{ 2 };
        CHECK(3 == *(s1 + s2));
        CHECK(-1 == *(s1 - s2));
        CHECK(2 == *(s1 * s2));
        CHECK(.5 == *(s1 / s2));
    }
}

TEST_CASE("Test LibPIC::CartCoord", "[LibPIC::CartCoord]")
{
    {
        constexpr CartCoord cart1{};
        CHECK(0 == cart1.x);

        constexpr CartCoord cart2{ 1 };
        CHECK(1 == cart2.x);

        constexpr auto cart3 = -cart2;
        CHECK(-1 == cart3.x);
    }
    {
        constexpr CartCoord cart1{ 1 };
        constexpr CartCoord cart2{ 2 };
        CHECK(1 == (+cart1).x);
        CHECK(-2 == (-cart2).x);
        CHECK(3 == (cart1 + cart2).x);
        CHECK(-1 == (cart1 - cart2).x);
        CHECK(2 == (cart1 * cart2).x);
        CHECK(.5 == (cart1 / cart2).x);
    }
}

TEST_CASE("Test LibPIC::CurviCoord", "[LibPIC::CurviCoord]")
{
    {
        constexpr CurviCoord curvi1{};
        CHECK(0 == curvi1.q1);

        constexpr CurviCoord curvi2{ 1 };
        CHECK(1 == curvi2.q1);

        constexpr auto curvi3 = -curvi2;
        CHECK(-1 == curvi3.q1);
    }
    {
        constexpr CurviCoord curvi1{ 1 };
        constexpr CurviCoord curvi2{ 2 };
        CHECK(1 == (+curvi1).q1);
        CHECK(-2 == (-curvi2).q1);
        CHECK(3 == (curvi1 + curvi2).q1);
        CHECK(-1 == (curvi1 - curvi2).q1);
        CHECK(2 == (curvi1 * curvi2).q1);
        CHECK(.5 == (curvi1 / curvi2).q1);
    }
}
