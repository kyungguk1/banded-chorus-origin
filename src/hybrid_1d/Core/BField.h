/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include "../ParamSet.h"

#include <HDF5Kit/HDF5Kit.h>

HYBRID1D_BEGIN_NAMESPACE
class EField;

class BField : public VectorGrid {
    VectorGrid B_prev;
    VectorGrid Ecovar;

public:
    ParamSet const params;
    Geometry const geomtr;

    explicit BField(ParamSet const &);
    BField &operator=(BField const &o) noexcept;

    void update(EField const &efield, Real dt) noexcept;

private:
    void mask(BField &B, MaskingFunction const &masking_function) const;
    void impl_update(BField &B_cart, VectorGrid const &Ecovar, Real cdtOsqrtg) const noexcept;
    auto cart_to_covar(VectorGrid &Ecovar, EField const &E_cart) const noexcept -> VectorGrid &;

    // attribute export facility
    //
    friend auto operator<<(hdf5::Group &obj, BField const &bfield) -> decltype(obj);
    friend auto operator<<(hdf5::Dataset &obj, BField const &bfield) -> decltype(obj);
    friend auto operator<<(hdf5::Group &&obj, BField const &bfield) -> decltype(obj)
    {
        return std::move(obj << bfield);
    }
    friend auto operator<<(hdf5::Dataset &&obj, BField const &bfield) -> decltype(obj)
    {
        return std::move(obj << bfield);
    }
};
HYBRID1D_END_NAMESPACE
