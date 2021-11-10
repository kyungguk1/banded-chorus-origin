/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "catch2/catch.hpp"

#include <Core/ExternalSource.h>
#include <exception>
#include <vector>

using namespace P1D;

TEST_CASE("Test pic_1d::ExternalSource", "[pic_1d::ExternalSource]")
{
    try {
        unsigned const rank   = 1;
        auto const     params = ParamSet{ rank, {} };

        Real const                  omega    = 2 * M_PI / (params.dt * 20);
        Real const                  start    = params.dt * 10;
        Real const                  duration = 2 * M_PI / omega * 5;
        Real const                  ease_in  = 2 * M_PI / omega * 3;
        unsigned const              N        = 3;
        ExternalSourceDesc<N> const desc{
            { omega, { start, duration }, ease_in },
            { ComplexVector{ 0, 1, 0 }, { 0, 0, -1i }, {} },
            { CurviCoord{ -3 }, CurviCoord{ 4 }, CurviCoord{ params.Nx } }
        };
        auto src = ExternalSource{ params, desc };
        src.set_cur_step(-1);
        REQUIRE(src.cur_step() == -1);

        // check envelope
        // 1. before ease-in
        for (long i = -100; i < -50; ++i) {
            REQUIRE(src.envelope(i * params.dt) == Approx{ 0 }.margin(1e-30));
        }
        // 2. ease-in phase
        for (long i = -50; i < 10; ++i) {
            auto const exact = .5 * (1 - std::cos(M_PI / 60 * Real(i + 50)));
            REQUIRE(src.envelope(i * params.dt) == Approx{ exact }.epsilon(1e-10));
        }
        // 3. middle phase
        for (long i = 10; i < 110; ++i) {
            auto const exact = 1;
            REQUIRE(src.envelope(i * params.dt) == Approx{ exact }.epsilon(1e-10));
        }
        // 4. ease-out phase
        for (long i = 110; i < 170; ++i) {
            auto const exact = .5 * (1 + std::cos(M_PI / 60 * Real(i - 110)));
            REQUIRE(src.envelope(i * params.dt) == Approx{ exact }.epsilon(1e-10));
        }
        // 5. before ease-out
        for (long i = 170; i < 200; ++i) {
            REQUIRE(src.envelope(i * params.dt) == Approx{ 0 }.margin(1e-30));
        }

        // check current
        Vector const J0re{ 1, M_SQRT1_2, 0 };
        Vector const J0im{ 0, M_SQRT1_2, 1 };
        for (long i = 10; i <= 30; ++i) {
            auto const t  = i * params.dt;
            auto const J  = src.current(J0re, J0im, t);
            auto const Jx = std::cos(2 * M_PI / 20 * Real(i - 10));
            REQUIRE(J.x == Approx{ Jx }.margin(1e-10));
            auto const Jz = std::sin(2 * M_PI / 20 * Real(i - 10));
            REQUIRE(J.z == Approx{ Jz }.margin(1e-10));
            auto const Jy = M_SQRT1_2 * (Jx + Jz);
            REQUIRE(J.y == Approx{ Jy }.margin(1e-10));
        }
    } catch (std::exception const &e) {
        INFO(e.what());
        CHECK(false);
    }
}
