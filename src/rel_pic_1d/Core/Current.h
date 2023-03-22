/*
 * Copyright (c) 2019-2022, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include "../ParamSet.h"

PIC1D_BEGIN_NAMESPACE
class Species;

/// current density
///
class Current : public Grid<CartVector> {
    Grid<CartVector> buffer;

public:
    ParamSet const params;

    explicit Current(ParamSet const &);
    Current &operator=(ParamSet const &) = delete;

    [[nodiscard]] auto &grid_whole_domain_extent() const noexcept { return params.half_grid_whole_domain_extent; }
    [[nodiscard]] auto &grid_subdomain_extent() const noexcept { return params.half_grid_subdomain_extent; }

    void reset() &noexcept { this->fill_all(CartVector{}); }
    void smooth() &noexcept { this->swap(buffer.smooth_assign(*this)); }

    Current &operator+=(Species const &) noexcept;
};
PIC1D_END_NAMESPACE
