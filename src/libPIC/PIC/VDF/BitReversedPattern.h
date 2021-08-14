/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC-config.h>

#include <limits>

LIBPIC_BEGIN_NAMESPACE
/// @brief Bit reversed pattern from Birdsall and Langdon (1985).
/// @discussion The original implementation is found in Kaijun Liu's PIC code.
///
/// The numbers will repeat once the `sequence` variable wraps around.
/// @note It satisfies the UniformRandomBitGenerator requirement.
///
template <unsigned base> class BitReversedPattern final {
    [[nodiscard]] static constexpr bool is_prime(unsigned const prime)
    {
        if (prime < 2)
            throw prime;
        unsigned i = prime;
        while (prime % --i) {}
        return 1 == i;
    }
    static_assert(base > 1 && is_prime(base), "base should be a prime number greater than 1");

public: // UniformRandomBitGenerator requirement
    using result_type = unsigned long;

    [[nodiscard]] static constexpr result_type min() noexcept { return 0; }
    [[nodiscard]] static constexpr result_type max() noexcept { return m_max; }

    [[nodiscard]] constexpr result_type operator()() noexcept { return next_pattern(m_seq++); }

public:
    BitReversedPattern(BitReversedPattern const &) = delete;
    BitReversedPattern &operator=(BitReversedPattern const &) = delete;
    constexpr BitReversedPattern() noexcept                   = default;

private:
    result_type                  m_seq{ 1 }; // sequence
    static constexpr result_type m_max = [x = result_type{ base }]() mutable noexcept {
        constexpr result_type max = std::numeric_limits<result_type>::max() / base;
        while (x < max) {
            x *= base;
        }
        return x; // base^n where n is an integer such that
                  // x < std::numeric_limits<result_type>::max()
    }();

    [[nodiscard]] static constexpr result_type next_pattern(result_type m_seq) noexcept
    {
        result_type power = max(), bit_pattern = 0;
        while (m_seq > 0) {
            bit_pattern += (m_seq % base) * (power /= base);
            m_seq /= base;
        }
        return bit_pattern;
    }
};
LIBPIC_END_NAMESPACE
