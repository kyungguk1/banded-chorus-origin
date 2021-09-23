/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "catch2/catch.hpp"

#include <PIC/CartCoord.h>
#include <PIC/CurviCoord.h>
#include <PIC/Scalar.h>

TEST_CASE("Test libPIC::Scalar", "[libPIC::Scalar]")
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

TEST_CASE("Test libPIC::CartCoord", "[libPIC::CartCoord]")
{
    {
        constexpr CartCoord s1{};
        CHECK(0 == *s1);
        CHECK(*s1 == s1());
        CHECK(*s1 == static_cast<double>(s1));

        constexpr CartCoord s2{ 1 };
        CHECK(1 == *s2);
        CHECK(*s2 == s2());
        CHECK(*s2 == static_cast<double>(s2));

        constexpr auto s3 = -s2;
        CHECK(-1 == *s3);
        CHECK(*s3 == s3());
        CHECK(*s3 == static_cast<double>(s3));
    }

    {
        constexpr CartCoord s1{ 1 };
        constexpr long      s2{ 2 };
        CHECK(3 == *(s1 + s2));
        CHECK(-1 == *(s1 - s2));
        CHECK(2 == *(s1 * s2));
        CHECK(.5 == *(s1 / s2));
    }
    {
        constexpr long      s1{ 1 };
        constexpr CartCoord s2{ 2 };
        CHECK(3 == *(s1 + s2));
        CHECK(-1 == *(s1 - s2));
        CHECK(2 == *(s1 * s2));
        CHECK(.5 == *(s1 / s2));
    }
}

TEST_CASE("Test libPIC::CurviCoord", "[libPIC::CurviCoord]")
{
    {
        constexpr CurviCoord s1{};
        CHECK(0 == *s1);
        CHECK(*s1 == s1());
        CHECK(*s1 == static_cast<double>(s1));

        constexpr CurviCoord s2{ 1 };
        CHECK(1 == *s2);
        CHECK(*s2 == s2());
        CHECK(*s2 == static_cast<double>(s2));

        constexpr auto s3 = -s2;
        CHECK(-1 == *s3);
        CHECK(*s3 == s3());
        CHECK(*s3 == static_cast<double>(s3));
    }

    {
        constexpr CurviCoord s1{ 1 };
        constexpr long       s2{ 2 };
        CHECK(3 == *(s1 + s2));
        CHECK(-1 == *(s1 - s2));
        CHECK(2 == *(s1 * s2));
        CHECK(.5 == *(s1 / s2));
    }
    {
        constexpr long       s1{ 1 };
        constexpr CurviCoord s2{ 2 };
        CHECK(3 == *(s1 + s2));
        CHECK(-1 == *(s1 - s2));
        CHECK(2 == *(s1 * s2));
        CHECK(.5 == *(s1 / s2));
    }
}
