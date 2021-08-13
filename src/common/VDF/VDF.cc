/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "VDF.h"

#include <algorithm>

COMMON_BEGIN_NAMESPACE
auto VDF::emit(unsigned int n) const -> std::vector<Particle>
{
    std::vector<Particle> ptls(n);
    std::generate(begin(ptls), end(ptls), [this]() {
        return this->emit();
    });
    return ptls;
}
COMMON_END_NAMESPACE
