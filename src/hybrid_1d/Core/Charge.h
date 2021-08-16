/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include "../ParamSet.h"

HYBRID1D_BEGIN_NAMESPACE
class Species;

/// charge density
///
class Charge : public ScalarGrid {
    ScalarGrid tmp;

public:
    ParamSet const params;

    virtual ~Charge() = default;
    explicit Charge(ParamSet const &);

    void reset() noexcept { this->fill(Scalar{ 0 }); }
    void smooth() noexcept
    {
        Grid::smooth(tmp, *this);
        this->swap(tmp);
    }

    virtual Charge &operator+=(Species const &sp) noexcept;
};

/// Λ
///
class Lambda : public Charge {
    using Charge::smooth;

public:
    using Charge::Charge;
    Lambda &operator+=(Species const &sp) noexcept override;
};
HYBRID1D_END_NAMESPACE
