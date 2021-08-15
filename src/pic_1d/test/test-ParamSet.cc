/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "catch2/catch.hpp"

#include <ParamSet.h>
#include <exception>

using namespace P1D;

TEST_CASE("Test pic_1d::ParamSet", "[pic_1d::ParamSet]")
{
    try {
        auto opts = Options{};
        opts.parse({ "--wd", "~/Downloads ", "-save ", "--load=false", "--outer_Nt=10" });

        auto const params = ParamSet{ 0, opts };
        CHECK(params.outer_Nt == 10);
        CHECK(params.working_directory == "~/Downloads");
        CHECK(params.snapshot_save == true);
        CHECK(params.snapshot_load == false);
        CHECK(params.domain_extent.loc == 0);
        CHECK(params.domain_extent.len == params.Nx / params.number_of_subdomains);
    } catch (std::exception const &e) {
        INFO(e.what());
        CHECK(false);
    }
}
