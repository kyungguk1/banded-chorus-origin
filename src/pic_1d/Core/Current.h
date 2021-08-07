/*
 * Copyright (c) 2019, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#ifndef Current_h
#define Current_h

#include "../Geometry.h"
#include "../ParamSet.h"

PIC1D_BEGIN_NAMESPACE
class Species;

/// current density
///
class Current : public VectorGrid {
    VectorGrid tmp;

public:
    ParamSet const params;
    Geometry const geomtr;

public:
    explicit Current(ParamSet const &);

    void reset() noexcept { this->fill(Vector{ 0 }); }
    void smooth() noexcept { _smooth(tmp, *this), this->swap(tmp); }

    Current &operator+=(Species const &sp) noexcept;
};
PIC1D_END_NAMESPACE

#endif /* Current_h */
