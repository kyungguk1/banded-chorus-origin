/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include "../ParamSet.h"

#include <HDF5Kit/HDF5Kit.h>

PIC1D_BEGIN_NAMESPACE
class BField;
class Current;

class EField : public VectorGrid {
    VectorGrid E_prev;
    VectorGrid Bcovar;

public:
    ParamSet const params;
    Geometry const geomtr;

    explicit EField(ParamSet const &);
    EField &operator=(EField const &) = delete;

    void update(BField const &bfield, Current const &current, Real dt) noexcept;

private:
    void mask(EField &E, MaskingFunction const &masking_function) const;
    void impl_update(EField &E_cart, VectorGrid const &B_covar, Real cdtOsqrtg, Current const &J_cart, Real dt) const noexcept;
    auto cart_to_covar(VectorGrid &Bcovar, BField const &B_cart) const noexcept -> VectorGrid &;

    // attribute export facility
    //
    friend auto operator<<(hdf5::Group &obj, EField const &efield) -> decltype(obj);
    friend auto operator<<(hdf5::Dataset &obj, EField const &efield) -> decltype(obj);
    friend auto operator<<(hdf5::Group &&obj, EField const &efield) -> decltype(obj)
    {
        return std::move(obj << efield);
    }
    friend auto operator<<(hdf5::Dataset &&obj, EField const &efield) -> decltype(obj)
    {
        return std::move(obj << efield);
    }
};
PIC1D_END_NAMESPACE
