/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC/Config.h>

#include <stdexcept>

LIBPIC_BEGIN_NAMESPACE
struct MaskingFunction {
    using Real = double;

    long masking_inset{}; // the length of the masking region near the one-side physical boundary
    Real masking_rate{};

    constexpr MaskingFunction() noexcept = default;
    constexpr MaskingFunction(unsigned masking_inset, Real masking_rate)
    : masking_inset{ masking_inset }, masking_rate{ masking_rate }
    {
        if (masking_rate < 0 || masking_rate > 1)
            throw std::invalid_argument{ __PRETTY_FUNCTION__ };
    }

    [[nodiscard]] constexpr Real operator()(Real const offset) const noexcept
    {
        if (0 == masking_inset || abs(offset) > masking_inset)
            return 1;

        auto const tmp = masking_rate * (1 - abs(offset) / masking_inset);
        return (1 - tmp) * (1 + tmp);
    }

private:
    [[nodiscard]] static constexpr Real abs(Real const x) noexcept { return x < 0 ? -x : x; }
};
LIBPIC_END_NAMESPACE
