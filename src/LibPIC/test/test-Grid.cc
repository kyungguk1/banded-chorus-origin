/*
 * Copyright (c) 2021-2022, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "catch2/catch.hpp"

#include <PIC/Grid.h>
#include <algorithm>
#include <numeric>

constexpr long Size = 10;

TEST_CASE("Test LibPIC::Grid_0", "[LibPIC::Grid_0]")
{
    constexpr long Pad = 0;
    using Grid         = PIC::Grid<Real, Size, Pad>;
    REQUIRE((Grid::size() == Size && Grid::pad_size() == Pad
             && Grid::max_size() == Grid::size() + Grid::pad_size() * 2));
    {
        // iterators and element access
        Grid g, &cg = g;
        REQUIRE(g.dead_begin() != nullptr);
        REQUIRE(std::distance(g.dead_begin(), g.dead_end()) == g.max_size());
        REQUIRE(std::distance(g.dead_begin(), g.begin()) == g.pad_size());
        REQUIRE(std::distance(g.end(), g.dead_end()) == g.pad_size());
        REQUIRE(std::distance(g.begin(), g.end()) == g.size());
        REQUIRE((g.dead_begin() == cg.dead_begin() && g.dead_end() == cg.dead_end()));
        REQUIRE((g.begin() == cg.begin() && g.end() == cg.end()));
        for (long i = -g.pad_size(); i < g.size() + g.pad_size(); ++i) {
            REQUIRE(&g[i] == std::next(g.begin(), i));
            REQUIRE(&cg[i] == &g[i]);
        }

        // fill
        CHECK(std::accumulate(g.dead_begin(), g.dead_end(), true, [](bool lhs, auto rhs) {
            return lhs && rhs == 0;
        }));
        g.fill(10);
        CHECK(std::accumulate(g.dead_begin(), g.dead_end(), true, [](bool lhs, auto rhs) {
            return lhs && rhs == 10;
        }));

        // copy/move/swap
        std::generate(g.dead_begin(), g.dead_end(), [i = 10]() mutable {
            return i++;
        });
        Grid g2;
        g2 = cg;
        REQUIRE(std::equal(cg.dead_begin(), cg.dead_end(), g2.dead_begin(), g2.dead_end()));
        g2.fill(10);
        g = std::move(g2);
        CHECK(std::accumulate(g.dead_begin(), g.dead_end(), true, [](bool lhs, auto rhs) {
            return lhs && rhs == 10;
        }));

        Grid g3{ std::move(g) };
        CHECK(std::accumulate(g3.dead_begin(), g3.dead_end(), true, [](bool lhs, auto rhs) {
            return lhs && rhs == 10;
        }));

        Grid g4;
        std::swap(g3, g4);
        CHECK(std::accumulate(g4.dead_begin(), g4.dead_end(), true, [](bool lhs, auto rhs) {
            return lhs && rhs == 10;
        }));
        CHECK(std::accumulate(g3.dead_begin(), g3.dead_end(), true, [](bool lhs, auto rhs) {
            return lhs && rhs == 0;
        }));
    }
}

TEST_CASE("Test LibPIC::Grid_1", "[LibPIC::Grid_1]")
{
    constexpr long Pad = 1;
    using Grid         = PIC::Grid<Real, Size, Pad>;
    REQUIRE((Grid::size() == Size && Grid::pad_size() == Pad
             && Grid::max_size() == Grid::size() + Grid::pad_size() * 2));
    {
        // iterators and element access
        Grid g, &cg = g;
        REQUIRE(g.dead_begin() != nullptr);
        REQUIRE(std::distance(g.dead_begin(), g.dead_end()) == g.max_size());
        REQUIRE(std::distance(g.dead_begin(), g.begin()) == g.pad_size());
        REQUIRE(std::distance(g.end(), g.dead_end()) == g.pad_size());
        REQUIRE(std::distance(g.begin(), g.end()) == g.size());
        REQUIRE((g.dead_begin() == cg.dead_begin() && g.dead_end() == cg.dead_end()));
        REQUIRE((g.begin() == cg.begin() && g.end() == cg.end()));
        for (long i = -g.pad_size(); i < g.size() + g.pad_size(); ++i) {
            REQUIRE(&g[i] == std::next(g.begin(), i));
            REQUIRE(&cg[i] == &g[i]);
        }

        // fill
        CHECK(std::accumulate(g.dead_begin(), g.dead_end(), true, [](bool lhs, auto rhs) {
            return lhs && rhs == 0;
        }));
        g.fill(10);
        CHECK(std::accumulate(g.dead_begin(), g.dead_end(), true, [](bool lhs, auto rhs) {
            return lhs && rhs == 10;
        }));

        // copy/move/swap
        std::generate(g.dead_begin(), g.dead_end(), [i = 10]() mutable {
            return i++;
        });
        Grid g2;
        g2 = cg;
        REQUIRE(std::equal(cg.dead_begin(), cg.dead_end(), g2.dead_begin(), g2.dead_end()));
        g2.fill(10);
        g = std::move(g2);
        CHECK(std::accumulate(g.dead_begin(), g.dead_end(), true, [](bool lhs, auto rhs) {
            return lhs && rhs == 10;
        }));

        Grid g3{ std::move(g) };
        CHECK(std::accumulate(g3.dead_begin(), g3.dead_end(), true, [](bool lhs, auto rhs) {
            return lhs && rhs == 10;
        }));

        Grid g4;
        std::swap(g3, g4);
        CHECK(std::accumulate(g4.dead_begin(), g4.dead_end(), true, [](bool lhs, auto rhs) {
            return lhs && rhs == 10;
        }));
        CHECK(std::accumulate(g3.dead_begin(), g3.dead_end(), true, [](bool lhs, auto rhs) {
            return lhs && rhs == 0;
        }));
    }
    { // smooth
        Grid filtered, source;
        source[4] = 1;
        CHECK(&filtered == &Grid::smooth(filtered, source));
        CHECK(1 == std::accumulate(filtered.dead_begin(), filtered.dead_end(), Real{}));
        std::accumulate(filtered.dead_begin(), std::next(filtered.begin() + 3), true,
                        [](bool lhs, auto rhs) {
                            return lhs && rhs == 0;
                        });
        std::accumulate(std::next(filtered.begin() + 6), filtered.dead_end(), true,
                        [](bool lhs, auto rhs) {
                            return lhs && rhs == 0;
                        });
        CHECK(.25 == filtered[3]);
        CHECK(.5 == filtered[4]);
        CHECK(.25 == filtered[5]);
    }
    { // Shape<1>
        using Shape = PIC::Shape<1>;
        Shape      sh;
        Real const weight = 10;
        Grid       g;

        g.fill(0);
        sh = Shape{ -1 };
        g.deposit(sh, weight);
        CHECK(std::accumulate(g.begin(), g.dead_end(), true, [](bool lhs, auto rhs) {
            return lhs && rhs == 0;
        }));
        CHECK(g[-1] == weight);

        g.fill(0);
        sh = Shape{ g.size() - 1e-15 }; // to avoid out-of-range memory access
        g.deposit(sh, weight);
        CHECK(std::accumulate(g.dead_begin(), g.end(), true, [](bool lhs, auto rhs) {
            return lhs && rhs == Approx{ 0 }.margin(1e-10);
        }));
        CHECK(*g.end() == Approx{ weight }.epsilon(1e-10));

        g.fill(0);
        sh = Shape{ 4.14 };
        g.deposit(sh, weight);
        CHECK(weight * sh.w<0>() == g[sh.i<0>()]);
        CHECK(weight * sh.w<1>() == g[sh.i<1>()]);
        CHECK(std::abs(weight - std::accumulate(g.dead_begin(), g.dead_end(), Real{})) / weight < 1e-15);

        g.fill(0);
        g[4] = weight;
        sh   = Shape{ 3.9 };
        CHECK(weight * sh.w<1>() == g.interp(sh));
        sh = Shape{ 4.1 };
        CHECK(weight * sh.w<0>() == g.interp(sh));
        sh = Shape{ 4.9 };
        CHECK(weight * sh.w<0>() == g.interp(sh));
    }
}

TEST_CASE("Test LibPIC::Grid_2", "[LibPIC::Grid_2]")
{
    constexpr long Pad = 2;
    using Grid         = PIC::Grid<Real, Size, Pad>;
    REQUIRE((Grid::size() == Size && Grid::pad_size() == Pad
             && Grid::max_size() == Grid::size() + Grid::pad_size() * 2));
    {
        // iterators and element access
        Grid g, &cg = g;
        REQUIRE(g.dead_begin() != nullptr);
        REQUIRE(std::distance(g.dead_begin(), g.dead_end()) == g.max_size());
        REQUIRE(std::distance(g.dead_begin(), g.begin()) == g.pad_size());
        REQUIRE(std::distance(g.end(), g.dead_end()) == g.pad_size());
        REQUIRE(std::distance(g.begin(), g.end()) == g.size());
        REQUIRE((g.dead_begin() == cg.dead_begin() && g.dead_end() == cg.dead_end()));
        REQUIRE((g.begin() == cg.begin() && g.end() == cg.end()));
        for (long i = -g.pad_size(); i < g.size() + g.pad_size(); ++i) {
            REQUIRE(&g[i] == std::next(g.begin(), i));
            REQUIRE(&cg[i] == &g[i]);
        }

        // fill
        CHECK(std::accumulate(g.dead_begin(), g.dead_end(), true, [](bool lhs, auto rhs) {
            return lhs && rhs == 0;
        }));
        g.fill(10);
        CHECK(std::accumulate(g.dead_begin(), g.dead_end(), true, [](bool lhs, auto rhs) {
            return lhs && rhs == 10;
        }));

        // copy/move/swap
        std::generate(g.dead_begin(), g.dead_end(), [i = 10]() mutable {
            return i++;
        });
        Grid g2;
        g2 = cg;
        REQUIRE(std::equal(cg.dead_begin(), cg.dead_end(), g2.dead_begin(), g2.dead_end()));
        g2.fill(10);
        g = std::move(g2);
        CHECK(std::accumulate(g.dead_begin(), g.dead_end(), true, [](bool lhs, auto rhs) {
            return lhs && rhs == 10;
        }));

        Grid g3{ std::move(g) };
        CHECK(std::accumulate(g3.dead_begin(), g3.dead_end(), true, [](bool lhs, auto rhs) {
            return lhs && rhs == 10;
        }));

        Grid g4;
        std::swap(g3, g4);
        CHECK(std::accumulate(g4.dead_begin(), g4.dead_end(), true, [](bool lhs, auto rhs) {
            return lhs && rhs == 10;
        }));
        CHECK(std::accumulate(g3.dead_begin(), g3.dead_end(), true, [](bool lhs, auto rhs) {
            return lhs && rhs == 0;
        }));
    }
    { // smooth
        Grid filtered, source;
        source[4] = 1;
        CHECK(&filtered == &Grid::smooth(filtered, source));
        CHECK(1 == std::accumulate(filtered.dead_begin(), filtered.dead_end(), Real{}));
        std::accumulate(filtered.dead_begin(), std::next(filtered.begin() + 3), true,
                        [](bool lhs, auto rhs) {
                            return lhs && rhs == 0;
                        });
        std::accumulate(std::next(filtered.begin() + 6), filtered.dead_end(), true,
                        [](bool lhs, auto rhs) {
                            return lhs && rhs == 0;
                        });
        CHECK(.25 == filtered[3]);
        CHECK(.5 == filtered[4]);
        CHECK(.25 == filtered[5]);
    }
    { // Shape<1>
        using Shape = PIC::Shape<1>;
        Shape      sh;
        Real const weight = 10;
        Grid       g;

        g.fill(0);
        sh = Shape{ -1 };
        g.deposit(sh, weight);
        CHECK(std::accumulate(g.begin(), g.dead_end(), true, [](bool lhs, auto rhs) {
            return lhs && rhs == 0;
        }));
        CHECK(g[-1] == weight);

        g.fill(0);
        sh = Shape{ g.size() };
        g.deposit(sh, weight);
        CHECK(std::accumulate(g.dead_begin(), g.end(), true, [](bool lhs, auto rhs) {
            return lhs && rhs == 0;
        }));
        CHECK(*g.end() == weight);

        g.fill(0);
        sh = Shape{ 4.14 };
        g.deposit(sh, weight);
        CHECK(weight * sh.w<0>() == g[sh.i<0>()]);
        CHECK(weight * sh.w<1>() == g[sh.i<1>()]);
        CHECK(std::abs(weight - std::accumulate(g.dead_begin(), g.dead_end(), Real{})) / weight
              < 1e-15);

        g.fill(0);
        g[4] = weight;
        sh   = Shape{ 3.9 };
        CHECK(weight * sh.w<1>() == g.interp(sh));
        sh = Shape{ 4.1 };
        CHECK(weight * sh.w<0>() == g.interp(sh));
        sh = Shape{ 4.9 };
        CHECK(weight * sh.w<0>() == g.interp(sh));
    }
    { // Shape<2>
        using Shape = PIC::Shape<2>;
        Shape      sh;
        Real const weight = 10;
        Grid       g;

        g.fill(0);
        sh = Shape{ -1 };
        g.deposit(sh, weight);
        CHECK(std::accumulate(std::next(g.begin()), g.dead_end(), true, [](bool lhs, auto rhs) {
            return lhs && rhs == 0;
        }));
        CHECK(g[-2] + g[-1] + g[0] == weight);

        g.fill(0);
        sh = Shape{ g.size() };
        g.deposit(sh, weight);
        CHECK(std::accumulate(g.dead_begin(), std::prev(g.end()), true, [](bool lhs, auto rhs) {
            return lhs && rhs == 0;
        }));
        CHECK(std::accumulate(std::prev(g.end()), g.dead_end(), Real{}) == weight);

        g.fill(0);
        sh = Shape{ 4.14 };
        g.deposit(sh, weight);
        CHECK(weight * sh.w<0>() == g[sh.i<0>()]);
        CHECK(weight * sh.w<1>() == g[sh.i<1>()]);
        CHECK(weight * sh.w<2>() == g[sh.i<2>()]);
        CHECK(std::abs(weight - std::accumulate(g.dead_begin(), g.dead_end(), Real{})) / weight
              < 1e-15);

        g.fill(0);
        g[4] = weight;
        sh   = Shape{ 3.9 };
        CHECK(weight * sh.w<1>() == g.interp(sh));
        sh = Shape{ 4.1 };
        CHECK(weight * sh.w<1>() == g.interp(sh));
        sh = Shape{ 4.9 };
        CHECK(weight * sh.w<0>() == g.interp(sh));
    }
}

TEST_CASE("Test LibPIC::Grid_3", "[LibPIC::Grid_3]")
{
    constexpr long Pad = 3;
    using Grid         = PIC::Grid<Real, Size, Pad>;
    REQUIRE((Grid::size() == Size && Grid::pad_size() == Pad
             && Grid::max_size() == Grid::size() + Grid::pad_size() * 2));
    {
        // iterators and element access
        Grid g, &cg = g;
        REQUIRE(g.dead_begin() != nullptr);
        REQUIRE(std::distance(g.dead_begin(), g.dead_end()) == g.max_size());
        REQUIRE(std::distance(g.dead_begin(), g.begin()) == g.pad_size());
        REQUIRE(std::distance(g.end(), g.dead_end()) == g.pad_size());
        REQUIRE(std::distance(g.begin(), g.end()) == g.size());
        REQUIRE((g.dead_begin() == cg.dead_begin() && g.dead_end() == cg.dead_end()));
        REQUIRE((g.begin() == cg.begin() && g.end() == cg.end()));
        for (long i = -g.pad_size(); i < g.size() + g.pad_size(); ++i) {
            REQUIRE(&g[i] == std::next(g.begin(), i));
            REQUIRE(&cg[i] == &g[i]);
        }

        // fill
        CHECK(std::accumulate(g.dead_begin(), g.dead_end(), true, [](bool lhs, auto rhs) {
            return lhs && rhs == 0;
        }));
        g.fill(10);
        CHECK(std::accumulate(g.dead_begin(), g.dead_end(), true, [](bool lhs, auto rhs) {
            return lhs && rhs == 10;
        }));

        // copy/move/swap
        std::generate(g.dead_begin(), g.dead_end(), [i = 10]() mutable {
            return i++;
        });
        Grid g2;
        g2 = cg;
        REQUIRE(std::equal(cg.dead_begin(), cg.dead_end(), g2.dead_begin(), g2.dead_end()));
        g2.fill(10);
        g = std::move(g2);
        CHECK(std::accumulate(g.dead_begin(), g.dead_end(), true, [](bool lhs, auto rhs) {
            return lhs && rhs == 10;
        }));

        Grid g3{ std::move(g) };
        CHECK(std::accumulate(g3.dead_begin(), g3.dead_end(), true, [](bool lhs, auto rhs) {
            return lhs && rhs == 10;
        }));

        Grid g4;
        std::swap(g3, g4);
        CHECK(std::accumulate(g4.dead_begin(), g4.dead_end(), true, [](bool lhs, auto rhs) {
            return lhs && rhs == 10;
        }));
        CHECK(std::accumulate(g3.dead_begin(), g3.dead_end(), true, [](bool lhs, auto rhs) {
            return lhs && rhs == 0;
        }));
    }
    { // smooth
        Grid filtered, source;
        source[4] = 1;
        CHECK(&filtered == &Grid::smooth(filtered, source));
        CHECK(1 == std::accumulate(filtered.dead_begin(), filtered.dead_end(), Real{}));
        std::accumulate(filtered.dead_begin(), std::next(filtered.begin() + 3), true,
                        [](bool lhs, auto rhs) {
                            return lhs && rhs == 0;
                        });
        std::accumulate(std::next(filtered.begin() + 6), filtered.dead_end(), true,
                        [](bool lhs, auto rhs) {
                            return lhs && rhs == 0;
                        });
        CHECK(.25 == filtered[3]);
        CHECK(.5 == filtered[4]);
        CHECK(.25 == filtered[5]);
    }
    { // Shape<1>
        using Shape = PIC::Shape<1>;
        Shape      sh;
        Real const weight = 10;
        Grid       g;

        g.fill(0);
        sh = Shape{ -1 };
        g.deposit(sh, weight);
        CHECK(std::accumulate(g.begin(), g.dead_end(), true, [](bool lhs, auto rhs) {
            return lhs && rhs == 0;
        }));
        CHECK(g[-1] == weight);

        g.fill(0);
        sh = Shape{ g.size() };
        g.deposit(sh, weight);
        CHECK(std::accumulate(g.dead_begin(), g.end(), true, [](bool lhs, auto rhs) {
            return lhs && rhs == 0;
        }));
        CHECK(*g.end() == weight);

        g.fill(0);
        sh = Shape{ 4.14 };
        g.deposit(sh, weight);
        CHECK(weight * sh.w<0>() == g[sh.i<0>()]);
        CHECK(weight * sh.w<1>() == g[sh.i<1>()]);
        CHECK(std::abs(weight - std::accumulate(g.dead_begin(), g.dead_end(), Real{})) / weight
              < 1e-15);

        g.fill(0);
        g[4] = weight;
        sh   = Shape{ 3.9 };
        CHECK(weight * sh.w<1>() == g.interp(sh));
        sh = Shape{ 4.1 };
        CHECK(weight * sh.w<0>() == g.interp(sh));
        sh = Shape{ 4.9 };
        CHECK(weight * sh.w<0>() == g.interp(sh));
    }
    { // Shape<2>
        using Shape = PIC::Shape<2>;
        Shape      sh;
        Real const weight = 10;
        Grid       g;

        g.fill(0);
        sh = Shape{ -1 };
        g.deposit(sh, weight);
        CHECK(std::accumulate(std::next(g.begin()), g.dead_end(), true, [](bool lhs, auto rhs) {
            return lhs && rhs == 0;
        }));
        CHECK(g[-2] + g[-1] + g[0] == weight);

        g.fill(0);
        sh = Shape{ g.size() };
        g.deposit(sh, weight);
        CHECK(std::accumulate(g.dead_begin(), std::prev(g.end()), true, [](bool lhs, auto rhs) {
            return lhs && rhs == 0;
        }));
        CHECK(std::accumulate(std::prev(g.end()), g.dead_end(), Real{}) == weight);

        g.fill(0);
        sh = Shape{ 4.14 };
        g.deposit(sh, weight);
        CHECK(weight * sh.w<0>() == g[sh.i<0>()]);
        CHECK(weight * sh.w<1>() == g[sh.i<1>()]);
        CHECK(weight * sh.w<2>() == g[sh.i<2>()]);
        CHECK(std::abs(weight - std::accumulate(g.dead_begin(), g.dead_end(), Real{})) / weight
              < 1e-15);

        g.fill(0);
        g[4] = weight;
        sh   = Shape{ 3.9 };
        CHECK(weight * sh.w<1>() == g.interp(sh));
        sh = Shape{ 4.1 };
        CHECK(weight * sh.w<1>() == g.interp(sh));
        sh = Shape{ 4.9 };
        CHECK(weight * sh.w<0>() == g.interp(sh));
    }
    { // Shape<3>
        using Shape = PIC::Shape<3>;
        Shape      sh;
        Real const weight = 10;
        Grid       g;

        g.fill(0);
        sh = Shape{ -1 };
        g.deposit(sh, weight);
        CHECK(std::accumulate(std::next(g.begin()), g.dead_end(), true, [](bool lhs, auto rhs) {
            return lhs && rhs == 0;
        }));
        CHECK(g[-3] + g[-2] + g[-1] + g[0] == weight);

        g.fill(0);
        sh = Shape{ g.size() };
        g.deposit(sh, weight);
        CHECK(std::accumulate(g.dead_begin(), std::prev(g.end(), 2), true, [](bool lhs, auto rhs) {
            return lhs && rhs == 0;
        }));
        CHECK(std::accumulate(std::prev(g.end(), 2), g.dead_end(), Real{}) == weight);

        g.fill(0);
        sh = Shape{ 4 + .1 };
        g.deposit(sh, weight);
        CHECK(weight * sh.w<0>() == g[sh.i<0>()]);
        CHECK(weight * sh.w<1>() == g[sh.i<1>()]);
        CHECK(weight * sh.w<2>() == g[sh.i<2>()]);
        CHECK(weight * sh.w<3>() == g[sh.i<3>()]);
        CHECK(std::abs(weight - std::accumulate(g.dead_begin(), g.dead_end(), Real{})) / weight
              < 1e-15);

        g.fill(0);
        g[4] = weight;
        sh   = Shape{ 3.9 };
        CHECK(weight * sh.w<2>() == g.interp(sh));
        sh = Shape{ 4.1 };
        CHECK(weight * sh.w<1>() == g.interp(sh));
        sh = Shape{ 4.9 };
        CHECK(weight * sh.w<1>() == g.interp(sh));
    }
}
