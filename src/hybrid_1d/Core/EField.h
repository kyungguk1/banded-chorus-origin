/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include "../ParamSet.h"

#include <HDF5Kit/HDF5Kit.h>

HYBRID1D_BEGIN_NAMESPACE
class BField;
class Charge;
class Current;

class EField : public VectorGrid {
    VectorGrid buffer;
    VectorGrid Je;
    VectorGrid dPe; // grad(Pe * c)
    //
    mutable ScalarGrid Pe;

public:
    ParamSet const params;
    Geometry const geomtr;

    explicit EField(ParamSet const &);
    EField &operator=(EField const &) = delete;

    void update(BField const &bfield, Charge const &charge, Current const &current) noexcept;

private:
    void mask(VectorGrid &grid, MaskingFunction const &masking_function) const;
    void impl_update_dPe(VectorGrid &grad_cPe_covar, Charge const &rho) const noexcept;
    void impl_update_Je(VectorGrid &Je_contr, Current const &Ji_cart, VectorGrid const &B_covar) const noexcept;
    void impl_update_E(EField &E_cart, Charge const &rho, VectorGrid const &dB_contr) const noexcept;

    auto cart_to_covar(VectorGrid &buffer, BField const &B_cart) const noexcept -> VectorGrid &;
    auto cart_to_contr(VectorGrid &buffer, BField const &B_cart) const noexcept -> VectorGrid &;

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
HYBRID1D_END_NAMESPACE
