/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "catch2/catch.hpp"

#include "../PIC/RandomReal.h"
#include "../PIC/splitmix64.h"
#include "../PIC/xoroshiro128.h"
#include <PIC/Predefined.h>
#include <algorithm>
#include <random>
#include <vector>

TEST_CASE("Test libPIC::splitmix64", "[libPIC::splitmix64]")
{
    constexpr auto seed     = 49847UL;
    constexpr auto a_number = splitmix64<seed>{}();
    CHECK(a_number > 0);

    auto rng     = splitmix64<seed>{};
    auto uniform = std::uniform_real_distribution<>{ 1, 10 };

    auto const        n_samples = 100000U;
    std::vector<Real> samples(n_samples);
    std::generate(begin(samples), end(samples), [&]() {
        return uniform(rng);
    });

    std::sort(begin(samples), end(samples));
    CHECK(samples.front() >= uniform.min());
    CHECK(samples.back() <= uniform.max());

    auto mean_sample = Real{};
    auto var_sample  = Real{};
    for (auto const &sample : samples) {
        mean_sample += sample / n_samples;
        var_sample += sample * sample / n_samples;
    }

    auto const mean_exact = (uniform.min() + uniform.max()) / 2;
    CHECK(std::abs(mean_sample - mean_exact) < mean_exact * 1e-2);
    auto const var_exact = (uniform.min() * uniform.min() + uniform.min() * uniform.max()
                            + uniform.max() * uniform.max())
                         / 3;
    CHECK(std::abs(var_sample - var_exact) < var_exact * 1e-2);
}

TEST_CASE("Test libPIC::xoroshiro128", "[libPIC::xoroshiro128]")
{

    constexpr auto seed     = 49847UL;
    constexpr auto a_number = xoroshiro128<seed>{}();
    CHECK(a_number > 0);

    auto rng     = xoroshiro128<seed>{};
    auto uniform = std::uniform_real_distribution<>{ 1, 10 };

    auto const        n_samples = 100000U;
    std::vector<Real> samples(n_samples);
    std::generate(begin(samples), end(samples), [&]() {
        return uniform(rng);
    });

    std::sort(begin(samples), end(samples));
    CHECK(samples.front() >= uniform.min());
    CHECK(samples.back() <= uniform.max());

    auto mean_sample = Real{};
    auto var_sample  = Real{};
    for (auto const &sample : samples) {
        mean_sample += sample / n_samples;
        var_sample += sample * sample / n_samples;
    }

    auto const mean_exact = (uniform.min() + uniform.max()) / 2;
    CHECK(std::abs(mean_sample - mean_exact) < mean_exact * 1e-3);
    auto const var_exact = (uniform.min() * uniform.min() + uniform.min() * uniform.max()
                            + uniform.max() * uniform.max())
                         / 3;
    CHECK(std::abs(var_sample - var_exact) < var_exact * 1e-3);
}

TEST_CASE("Test libPIC::RandomReal::uniform_real", "[libPIC::RandomReal::uniform_real]")
{

    constexpr auto seed = 49847UL;
    constexpr Real min = 1, max = 10;

    auto const        n_samples = 100000U;
    std::vector<Real> samples(n_samples);
    std::generate(begin(samples), end(samples), []() {
        return uniform_real<seed>() * (max - min) + min;
    });

    std::sort(begin(samples), end(samples));
    CHECK(samples.front() >= min);
    CHECK(samples.back() <= max);

    auto mean_sample = Real{};
    auto var_sample  = Real{};
    for (auto const &sample : samples) {
        mean_sample += sample / n_samples;
        var_sample += sample * sample / n_samples;
    }

    auto const mean_exact = (min + max) / 2;
    CHECK(std::abs(mean_sample - mean_exact) < mean_exact * 1e-2);
    auto const var_exact = (min * min + min * max + max * max) / 3;
    CHECK(std::abs(var_sample - var_exact) < var_exact * 1e-2);
}

TEST_CASE("Test libPIC::RandomReal::bit_reversed", "[libPIC::RandomReal::bit_reversed]")
{

    constexpr auto base = 11U;
    constexpr Real min = 1, max = 10;

    auto const        n_samples = 100000U;
    std::vector<Real> samples(n_samples);
    std::generate(begin(samples), end(samples), []() {
        return bit_reversed<base>() * (max - min) + min;
    });

    std::sort(begin(samples), end(samples));
    CHECK(samples.front() >= min);
    CHECK(samples.back() <= max);

    auto mean_sample = Real{};
    auto var_sample  = Real{};
    for (auto const &sample : samples) {
        mean_sample += sample / n_samples;
        var_sample += sample * sample / n_samples;
    }

    auto const mean_exact = (min + max) / 2;
    CHECK(std::abs(mean_sample - mean_exact) < mean_exact * 1e-4);
    auto const var_exact = (min * min + min * max + max * max) / 3;
    CHECK(std::abs(var_sample - var_exact) < var_exact * 1.1e-4);
}
