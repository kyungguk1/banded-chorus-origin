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
        constexpr auto badge1() const noexcept { return Badge<S>{}; }
        // constexpr auto badge2() const noexcept { return Badge<long>{}; }
    };
    [[maybe_unused]] constexpr auto badge = S{}.badge1();
    REQUIRE_NOTHROW(S{}.badge1());
    // REQUIRE_NOTHROW(S{}.badge2());
}
