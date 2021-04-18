//
//  BitReversedPattern.cc
//  hybrid_1d
//
//  Created by KYUNGGUK MIN on 11/30/19.
//  Copyright © 2019 Kyungguk Min & Kaijun Liu. All rights reserved.
//

#include "BitReversedPattern.h"

#include "../Utility/println.h"

#include <iostream>

using H1D::BitReversedPattern;

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
void H1D::test_BitReversedPattern()
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
void H1D::test_BitReversedPattern()
{
}
#endif
