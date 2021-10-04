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

        REQUIRE_THROWS_AS(ParamSet(ParamSet::number_of_subdomains, opts), std::invalid_argument);

        unsigned const rank   = 1;
        auto const     params = ParamSet{ rank, opts };
        CHECK(params.outer_Nt == 10);
        CHECK(params.working_directory == "~/Downloads");
        CHECK(params.snapshot_save == true);
        CHECK(params.snapshot_load == false);

        CHECK(params.full_grid_whole_domain_extent.loc == -.5 * params.Nx);
        CHECK(params.full_grid_whole_domain_extent.len == params.Nx);
        CHECK(params.half_grid_whole_domain_extent.loc == params.full_grid_whole_domain_extent.min() + 0.5);
        CHECK(params.half_grid_whole_domain_extent.len == params.Nx);

        Real const Mx = params.Nx / params.number_of_subdomains;
        CHECK(params.full_grid_subdomain_extent.loc == params.full_grid_whole_domain_extent.min() + rank * Mx);
        CHECK(params.full_grid_subdomain_extent.len == Mx);
        CHECK(params.half_grid_subdomain_extent.loc == params.full_grid_subdomain_extent.min() + 0.5);
        CHECK(params.half_grid_subdomain_extent.len == Mx);
    } catch (std::exception const &e) {
        INFO(e.what());
        CHECK(false);
    }
}
