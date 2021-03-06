/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <ostream>
#include <utility>

// synthetic sugar for output stream operator (<<)
//
namespace {

template <class CharT, class Traits, class... Args>
decltype(auto) print(std::basic_ostream<CharT, Traits> &os, Args &&...args)
{
    static_assert(sizeof...(Args) > 0, "there should be at least one item to print");
    return (os << ... << args);
}
//
template <class CharT, class Traits, class... Args>
decltype(auto) println(std::basic_ostream<CharT, Traits> &os, Args &&...args)
{
    static_assert(sizeof...(Args) > 0, "there should be at least one item to print");
    return print(os, std::forward<Args>(args)...) << '\n';
}

} // namespace
