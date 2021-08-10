/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "catch2/catch.hpp"

#include <Utility/Options.h>
#include <Utility/println.h>
#include <algorithm>
#include <iostream>
#include <iterator>

using COMMON_NAMESPACE::Options;

TEST_CASE("Test common::Options", "[common::Options]")
{
    println(std::cout, "in ", __PRETTY_FUNCTION__);

    Options opts;
    opts.parse({ "--save=false", "--long=3", "--dir", "~" });
    auto const unparsed = opts.parse({ { "a", "- save  ", "b", "-", "--", "--load=false",
                                         "--long = -3", "--str= s", "-abc xyz" } });

    if (!unparsed.empty()) {
        print(std::cout, "unparsed arguments = ");
        print(std::cout, unparsed.front());
        std::for_each(std::next(begin(unparsed)), end(unparsed), [](auto const &arg) {
            print(std::cout, ", ", arg);
        });
        std::cout << '\n';
    }
    println(std::cout, "options = ", opts);
    println(std::cout, "str = ", opts->at("str").as<char const *>());
    println(std::cout, "save = ", opts->at("save").as<bool>());
    println(std::cout, "load = ", opts->at("load").as<bool>());
    println(std::cout, "long = ", opts->at("long").as<long>());
    println(std::cout, "dir = ", opts->at("dir").as<char const *>());
    println(std::cout, "abc xyz = ", opts->at("abc xyz").as<std::string>());

    std::vector<std::string> const unparsible{ "a", "b", "-", "--" };
    CHECK(std::equal(begin(unparsed), end(unparsed), begin(unparsible), end(unparsible)));

    std::map<std::string, std::string> parsible_opts{ { "abc xyz", "true" }, { "dir", "~" },
                                                      { "load", "false" },   { "long", "-3" },
                                                      { "save", "true" },    { "str", "s" } };
    REQUIRE(std::equal(begin(*opts), end(*opts), begin(parsible_opts), end(parsible_opts),
                       [](auto const &lhs, auto const &rhs) {
                           return lhs.first == rhs.first && *lhs.second == rhs.second;
                       }));
    REQUIRE(parsible_opts.at("str") == opts->at("str").as<char const *>());
    REQUIRE(opts->at("save").as<bool>());
    REQUIRE_FALSE(opts->at("load").as<bool>());
    REQUIRE(-3 == opts->at("long").as<long>());
}
