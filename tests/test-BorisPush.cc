/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "catch2/catch.hpp"

#include <Utility/BorisPush.h>
#include <Utility/println.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <numeric>
#include <vector>

using COMMON_NAMESPACE::BorisPush;
using COMMON_NAMESPACE::Vector;
using Real = COMMON_NAMESPACE::Real;

namespace {
[[maybe_unused]] void printer(std::vector<Vector> const &vs)
{
    if (vs.empty()) {
        std::puts("{}");
    } else {
        print(std::cout, "{\n    ", vs.front());
        std::for_each(std::next(begin(vs)), end(vs), [](Vector const &v) {
            print(std::cout, ",\n    ", v);
        });
        std::puts("\n}");
    }
}
[[maybe_unused]] void printer(std::vector<std::pair<Real, Vector>> const &gvs)
{
    if (gvs.empty()) {
        std::puts("{}");
    } else {
        print(std::cout, "{\n    {", gvs.front().first, ", ", gvs.front().second, '}');
        std::for_each(std::next(begin(gvs)), end(gvs), [](std::pair<Real, Vector> const &v) {
            print(std::cout, ",\n    {", v.first, ", ", v.second, '}');
        });
        std::puts("\n}");
    }
}
} // namespace

TEST_CASE("Test common::NonrelativisticBorisPush", "[common::NonrelativisticBorisPush]")
{
    unsigned const  nt = 360;
    Real const      dt = 1;
    Real const      c  = 1;
    Real const      O0 = 1;
    Real const      Oc = 2 * M_PI / nt;
    BorisPush const boris{ dt, c, O0, Oc };

    std::vector<Vector> vs;
    std::generate_n(std::back_inserter(vs), nt,
                    [&boris, v = Vector{ 1, 0, 0 }]() mutable -> Vector {
                        constexpr Vector B0{ 0, 0, 1 };
                        boris.non_relativistic(v, B0, {});
                        return v;
                    });
    // printer(vs);
    CHECK(std::accumulate(begin(vs), end(vs), true, [](bool lhs, Vector const &v) {
        return lhs && std::abs(dot(v, v) - 1) < 1e-15 && v.z == 0;
    }));

    vs = { { 0, 0, 0 } };
    for (long t = 0; t < nt; ++t) {
        Vector const E{ 0, 0, std::cos(2 * M_PI * (Real(t) + dt / 2) / nt) };
        boris.non_relativistic(vs.emplace_back(vs.back()), {}, E / c);
    }
    vs.erase(begin(vs));
    // printer(vs);
    CHECK(std::accumulate(begin(vs), end(vs), true, [](bool lhs, Vector const &v) {
        return lhs && v.x == 0 && v.y == 0;
    }));
    auto const rms = std::accumulate(begin(vs), end(vs), Real{},
                                     [](Real lhs, Vector const &v) {
                                         return lhs += dot(v, v);
                                     })
                   / Real(vs.size());
    CHECK(std::abs(0.50001269258580949284 - rms) / rms < 1e-10);
}

TEST_CASE("Test common::RelativisticBorisPush", "[common::RelativisticBorisPush]")
{
    { // non-relativistic regime
        unsigned const  nt = 360;
        Real const      dt = 1;
        Real const      c  = 1e3;
        Real const      O0 = 1;
        Real const      Oc = 2 * M_PI / nt;
        BorisPush const boris{ dt, c, O0, Oc };

        std::vector<std::pair<Real, Vector>> gvs;
        std::generate_n(std::back_inserter(gvs), nt, [&boris, gv = Vector{ 1, 0, 0 }]() mutable {
            constexpr Vector B0{ 0, 0, 1 };
            auto             g = boris.relativistic(gv, B0, {});
            return std::make_pair(g, gv / g);
        });
        // printer(gvs);
        CHECK(std::accumulate(begin(gvs), end(gvs), true,
                              [gamma = std::sqrt(1 + 1 / (c * c)), c](bool lhs, auto const &gv) {
                                  auto const &[_, v] = gv;
                                  auto beta          = std::sqrt(dot(v, v)) / c;
                                  return lhs
                                      && std::abs(1 / std::sqrt((1 - beta) * (1 + beta)) - gamma)
                                             < 1e-15
                                      && v.z == 0;
                              }));
        CHECK(std::accumulate(begin(gvs), end(gvs), true,
                              [gamma = std::sqrt(1 + 1 / (c * c))](bool lhs, auto const &gv) {
                                  auto const &[g, _] = gv;
                                  return lhs && std::abs(g - gamma) / gamma < 1e-15;
                              }));

        gvs = { { 0, { 0, 0, 0 } } };
        for (long t = 0; t < nt; ++t) {
            Vector const E{ 0, 0, std::cos(2 * M_PI * (Real(t) + dt / 2) / nt) };
            auto        &gv = gvs.emplace_back(gvs.back());
            gv.first        = boris.relativistic(gv.second, {}, E / c);
        }
        gvs.erase(begin(gvs));
        // printer(gvs);
        CHECK(std::accumulate(begin(gvs), end(gvs), true, [](bool lhs, auto const &gv) {
            auto const &[g, v] = gv;
            return lhs && v.x == 0 && v.y == 0;
        }));
        auto const rms = std::accumulate(begin(gvs), end(gvs), Real{},
                                         [](Real lhs, auto const &gv) {
                                             auto const &[g, v] = gv;
                                             return lhs += dot(v, v);
                                         })
                       / Real(gvs.size());
        CHECK(std::abs(0.50001269258580949284 - rms) / rms < 1e-10);
    }
    { // relativistic regime
        unsigned const  nt = 360;
        Real const      dt = 1;
        Real const      c  = 4;
        Real const      O0 = 1;
        Real const      Oc = 2 * M_PI / nt;
        BorisPush const boris{ dt, c, O0, Oc };

        std::vector<std::pair<Real, Vector>> gvs;
        std::generate_n(std::back_inserter(gvs), nt, [&boris, gv = Vector{ 1, 0, 0 }]() mutable {
            constexpr Vector B0{ 0, 0, 1 };
            auto             g = boris.relativistic(gv, B0, {});
            return std::make_pair(g, gv / g);
        });
        // printer(gvs);
        CHECK(std::accumulate(begin(gvs), end(gvs), true,
                              [gamma = std::sqrt(1 + 1 / (c * c)), c](bool lhs, auto const &gv) {
                                  auto const &[_, v] = gv;
                                  auto beta          = std::sqrt(dot(v, v)) / c;
                                  return lhs
                                      && std::abs(1 / std::sqrt((1 - beta) * (1 + beta)) - gamma)
                                             < 1e-15
                                      && v.z == 0;
                              }));
        CHECK(std::accumulate(begin(gvs), end(gvs), true,
                              [gamma = std::sqrt(1 + 1 / (c * c))](bool lhs, auto const &gv) {
                                  auto const &[g, _] = gv;
                                  return lhs && std::abs(g - gamma) / gamma < 1e-15;
                              }));

        gvs = { { 0, { 0, 0, 0 } } };
        for (long t = 0; t < nt; ++t) {
            Vector const E{ 0, 0, std::cos(2 * M_PI * (Real(t) + dt / 2) / nt) };
            auto        &gv = gvs.emplace_back(gvs.back());
            gv.first        = boris.relativistic(gv.second, {}, E / c);
        }
        gvs.erase(begin(gvs));
        // printer(gvs);
        CHECK(std::accumulate(begin(gvs), end(gvs), true, [](bool lhs, auto const &gv) {
            auto const &[g, v] = gv;
            return lhs && v.x == 0 && v.y == 0;
        }));
        auto const rms = std::accumulate(begin(gvs), end(gvs), Real{},
                                         [](Real lhs, auto const &gv) {
                                             auto const &[g, v] = gv;
                                             return lhs += dot(v, v);
                                         })
                       / Real(gvs.size());
        CHECK(std::abs(0.50001269258580949284 - rms) / rms < 1e-10);
    }
}
