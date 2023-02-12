/*
 * Copyright (c) 2023, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "test-Math.h"
#include "catch2/catch.hpp"
#include <PIC/Math/Faddeeva.hh>

TEST_CASE("Test LibPIC::Math::Dawson", "[LibPIC::Math::Dawson]")
{
    long i = 0;
    for (auto const &[x, f1] : dawson_F) {
        auto const f2 = Faddeeva::Dawson(x);
        INFO("i = " << i << ", x = " << x << ", f1 = " << f1 << ", f2 = " << f2);
        REQUIRE(f1 == Approx{ f2 }.epsilon(1e-10));
        ++i;
    }
}
