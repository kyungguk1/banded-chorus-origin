/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "catch2/catch.hpp"

#include <PIC/Badge.h>

TEST_CASE("Test libPIC::Badge", "[libPIC::Badge]")
{
    struct S {
        auto badge() const noexcept { return Badge<S>{}; }
    };
    REQUIRE_NOTHROW(S{}.badge());
}
