/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#ifndef EField_h
#define EField_h

#include "../Geometry.h"
#include "../ParamSet.h"

#include <HDF5Kit/HDF5Kit.h>

PIC1D_BEGIN_NAMESPACE
class BField;
class Current;

class EField : public VectorGrid {
public:
    ParamSet const params;
    Geometry const geomtr;

public:
    explicit EField(ParamSet const &);

    void update(BField const &bfield, Current const &current, Real dt) noexcept;

private:
    static inline void _update(EField &E, BField const &B, Real cdtODx, Current const &J,
                               Real dt) noexcept;

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

#endif /* EField_h */
