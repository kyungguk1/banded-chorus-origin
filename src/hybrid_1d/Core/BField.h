/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#ifndef BField_h
#define BField_h

#include "../Geometry.h"
#include "../ParamSet.h"

#include <HDF5Kit/HDF5Kit.h>

HYBRID1D_BEGIN_NAMESPACE
class EField;

class BField : public VectorGrid {
public:
    ParamSet const params;
    Geometry const geomtr;

public:
    explicit BField(ParamSet const &);
    BField &operator=(BField const &o) noexcept
    {
        this->GridQ::operator=(o);
        return *this;
    }
    using GridQ::swap;

    void update(EField const &efield, Real dt) noexcept;

private:
    static inline void _update(BField &B, EField const &E, Real cdtODx) noexcept;

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

#endif /* BField_h */
