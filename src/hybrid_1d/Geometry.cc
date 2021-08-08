/*
 * Copyright (c) 2019, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "Geometry.h"

#include <cmath>

H1D::Geometry::Geometry(Input const &params) noexcept
{
    Real const theta = params.theta * M_PI / 180;
    this->e1         = Vector{ std::cos(theta), std::sin(theta), 0 };
    this->e2         = Vector{ -std::sin(theta), std::cos(theta), 0 };
    this->B0         = this->e1 * params.O0;
}
