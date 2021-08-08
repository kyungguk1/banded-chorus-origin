/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "BField.h"

#include "./EField.h"

H1D::BField::BField(ParamSet const &params) : GridQ{}, params{ params }, geomtr{ params }
{
    this->fill(geomtr.B0); // fill with background B
}

void H1D::BField::update(EField const &efield, Real const dt) noexcept
{
    Real const cdtODx = dt * params.c / params.Dx;
    _update(*this, efield, cdtODx);
}
void H1D::BField::_update(BField &B, EField const &E, Real const cdtODx) noexcept
{
    for (long i = 0; i < E.size(); ++i) {
        B[i].x += 0;
        B[i].y += (+E[i - 0].z - E[i - 1].z) * cdtODx;
        B[i].z += (-E[i - 0].y + E[i - 1].y) * cdtODx;
    }
}

namespace {
template <class Object>
decltype(auto) write_attr(Object &obj, [[maybe_unused]] H1D::BField const &bfield)
{
    return obj;
}
} // namespace
auto H1D::operator<<(hdf5::Group &obj, [[maybe_unused]] H1D::BField const &bfield) -> decltype(obj)
{
    return write_attr(obj, bfield);
}
auto H1D::operator<<(hdf5::Dataset &obj, [[maybe_unused]] H1D::BField const &bfield)
    -> decltype(obj)
{
    return write_attr(obj, bfield);
}
