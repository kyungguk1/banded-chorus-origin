/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include "../ParamSet.h"

PIC1D_BEGIN_NAMESPACE
class Species;

/// current density
///
class Current : public VectorGrid {
    VectorGrid buffer;

public:
    ParamSet const params;

    explicit Current(ParamSet const &);
    Current &operator=(ParamSet const &) = delete;

    void reset() noexcept { this->fill(Vector{ 0 }); }
    void smooth() noexcept
    {
        Grid::smooth(buffer, *this);
        this->swap(buffer);
    }

    Current &operator+=(Species const &sp) noexcept;
};
PIC1D_END_NAMESPACE
