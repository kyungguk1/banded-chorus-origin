/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "catch2/catch.hpp"

#include <Utility/Scalar.h>

using COMMON_NAMESPACE::Scalar;

TEST_CASE("Test common::Scalar", "[common::Scalar]")
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
