/*
 * Copyright (c) 2019, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#ifndef Charge_h
#define Charge_h

#include "../Geometry.h"
#include "../ParamSet.h"

HYBRID1D_BEGIN_NAMESPACE
class Species;

/// charge density
///
class Charge : public ScalarGrid {
    ScalarGrid tmp;

public:
    ParamSet const params;
    Geometry const geomtr;

public:
    virtual ~Charge() = default;
    explicit Charge(ParamSet const &);

    void reset() noexcept { this->fill(Scalar{ 0 }); }
    void smooth() noexcept { _smooth(tmp, *this), this->swap(tmp); }

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

#endif /* Charge_h */
