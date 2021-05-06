/*
 * Copyright (c) 2019, Kyungguk Min
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * - Redistributions of source code must retain the above copyright notice, this
 *  list of conditions and the following disclaimer.
 *
 * - Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef BitReversedPattern_h
#define BitReversedPattern_h

#include "../Macros.h"

#include <limits>

PIC1D_BEGIN_NAMESPACE
/**
 @brief Bit reversed pattern from Birdsall and Langdon (1985).
 @discussion The original implementation is found in Kaijun's PIC code.

 The numbers will repeat once the `sequence' variable wraps around.
 @note It satisfies the UniformRandomBitGenerator requirement.
 */
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
    [[nodiscard]] static constexpr result_type max() noexcept { return _max; }

    [[nodiscard]] result_type operator()() noexcept { return next_pattern(sequence++); }

public:
    BitReversedPattern(BitReversedPattern const &) = delete;
    BitReversedPattern &operator=(BitReversedPattern const &) = delete;
    constexpr BitReversedPattern() noexcept                   = default;

private:
    result_type                  sequence{ 1 };
    static constexpr result_type _max = [x = result_type{ base }]() mutable noexcept {
        constexpr result_type max = std::numeric_limits<result_type>::max() / base;
        while (x < max) {
            x *= base;
        }
        return x; // base^n where n is an integer such that
                  // x < std::numeric_limits<result_type>::max()
    }();

    [[nodiscard]] static constexpr result_type next_pattern(result_type sequence) noexcept
    {
        result_type power = max(), bit_pattern = 0;
        while (sequence > 0) {
            bit_pattern += (sequence % base) * (power /= base);
            sequence /= base;
        }
        return bit_pattern;
    }
};

// not for public use
//
void test_BitReversedPattern();
PIC1D_END_NAMESPACE

#endif /* BitReversedPattern_h */
