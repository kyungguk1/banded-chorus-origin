/*
 * Copyright (c) 2020-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC/Config.h>
#include <PIC/println.h>

#include <sstream>
#include <string>
#include <utility>

// error handling routines
//
LIBPIC_BEGIN_NAMESPACE
void print_backtrace();

[[noreturn, maybe_unused]] void fatal_error(char const *reason) noexcept;
[[noreturn, maybe_unused]] void fatal_error(std::string const &reason) noexcept;
[[noreturn, maybe_unused]] void lippincott() noexcept;

template <class... Args>
[[noreturn, maybe_unused]] void fatal_error(Args &&...args)
{
    std::ostringstream ss;
    {
        print(ss, std::forward<Args>(args)...);
    }
    std::string const str = ss.str(); // should hold onto it
    fatal_error(str);
}
LIBPIC_END_NAMESPACE
