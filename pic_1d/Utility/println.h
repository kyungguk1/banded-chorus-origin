//
//  println.h
//  pic_1d
//
//  Created by KYUNGGUK MIN on 1/30/19.
//  Copyright © 2019 Kyungguk Min & Kaijun Liu. All rights reserved.
//

#ifndef println_h
#define println_h

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

#endif /* println_h */
