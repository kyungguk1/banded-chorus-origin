/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include "../ParamSet.h"

HYBRID1D_BEGIN_NAMESPACE
class BField;
class EField;
class Lambda;
class Species;
class Gamma;

/// current density
///
class Current : public VectorGrid {
    VectorGrid buffer;

public:
    ParamSet const params;
    Geometry const geomtr;

    virtual ~Current() = default;
    explicit Current(ParamSet const &);
    Current &operator=(ParamSet const &) = delete;

    void reset() noexcept { this->fill(Vector{ 0 }); }
    void smooth() noexcept
    {
        Grid::smooth(buffer, *this);
        this->swap(buffer);
    }

    virtual Current &operator+=(Species const &sp) noexcept;

    void advance(Lambda const &lambda, Gamma const &gamma, BField const &bfield, EField const &efield, Real dt) noexcept;

private:
    void impl_advance(Current &J, Lambda const &L, Gamma const &G, BField const &dB, EField const &E, Real dt) const noexcept;
};

/// Γ
///
class Gamma : public Current {
    using Current::advance;
    using Current::smooth;

public:
    using Current::Current;
    Gamma &operator+=(Species const &sp) noexcept override;
};
HYBRID1D_END_NAMESPACE
