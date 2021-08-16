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
    VectorGrid Je;
    ScalarGrid Pe;

public:
    ParamSet const params;

    explicit EField(ParamSet const &);

    void update(BField const &bfield, Charge const &charge, Current const &current) noexcept;

private:
    void impl_update_Pe(ScalarGrid &Pe, Charge const &rho) const noexcept;
    void impl_update_Je(VectorGrid &Je, Current const &Ji, BField const &B) const noexcept;
    void impl_update_E(EField &E, BField const &B, Charge const &rho) const noexcept;

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
