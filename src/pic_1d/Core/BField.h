/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include "../ParamSet.h"
#include <PIC/Geometry.h>

#include <HDF5Kit/HDF5Kit.h>

PIC1D_BEGIN_NAMESPACE
class EField;

class BField : public VectorGrid {
public:
    ParamSet const params;
    Geometry const geomtr;

    explicit BField(ParamSet const &);
    BField &operator=(BField const &o) noexcept
    { // this is only for the return type casting
        this->Grid::operator=(o);
        return *this;
    }

    void update(EField const &efield, Real dt) noexcept;

private:
    static void impl_update(BField &B, EField const &E, Real cdtODx) noexcept;

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
PIC1D_END_NAMESPACE
