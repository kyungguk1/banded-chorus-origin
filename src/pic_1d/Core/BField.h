/*
 * Copyright (c) 2019-2021, Kyungguk Min
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * - Redistributions of source code must retain the above copyright notice, this
 *  list of conditions and the following disclaimer.
 *
 * - Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef BField_h
#define BField_h

#include "../Geometry.h"
#include "../ParamSet.h"

#include <HDF5Kit/HDF5Kit.h>

PIC1D_BEGIN_NAMESPACE
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
PIC1D_END_NAMESPACE

#endif /* BField_h */
