/*
 * Copyright (c) 2020-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "lippincott.h"

#include <algorithm>
#include <array>
#include <cstdio>
#include <execinfo.h>

LIBPIC_BEGIN_NAMESPACE
void print_backtrace()
{
    constexpr unsigned             stack_size = 64;
    std::array<void *, stack_size> array{};
    int const                      size = backtrace(array.data(), std::size(array));
    if (char **const strings = backtrace_symbols(array.data(), size)) {
        std::for_each_n(strings, size, &std::puts);
        free(strings);
    }
}
LIBPIC_END_NAMESPACE
