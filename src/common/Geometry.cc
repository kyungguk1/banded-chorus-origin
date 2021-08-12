/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "Geometry.h"

#include <cmath>
#include <stdexcept>

COMMON_BEGIN_NAMESPACE
Geometry::Geometry(const Vector &_B0) : B0{ _B0 }
{
    auto const mag = std::sqrt(dot(B0, B0));
    if (std::abs(B0.z) > mag * 1e-15)
        throw std::invalid_argument{ __PRETTY_FUNCTION__ };

    B0.z = 0;
    e1   = B0 / mag;
    e2   = cross(e3, e1);
}
COMMON_END_NAMESPACE
