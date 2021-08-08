/*
 * Copyright (c) 2019, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "BitReversedPattern.h"

#include "../Utility/println.h"

#include <iostream>

using P1D::BitReversedPattern;

#if defined(DEBUG)
namespace {
template <unsigned base> void test(BitReversedPattern<base> g)
{
    print(std::cout, '{', g());
    for (long i = 1; i < 200; ++i) {
        print(std::cout, ", ", g());
    }
    println(std::cout, '}');
}
} // namespace
void P1D::test_BitReversedPattern()
{
    print(std::cout, '{');
    test<2>({});
    print(std::cout, ", ");
    test<3>({});
    print(std::cout, ", ");
    test<5>({});
    print(std::cout, ", ");
    test<19>({});
    println(std::cout, '}');
}
#else
void P1D::test_BitReversedPattern()
{
}
#endif
