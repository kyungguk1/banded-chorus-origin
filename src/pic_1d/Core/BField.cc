/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "BField.h"
#include "EField.h"

PIC1D_BEGIN_NAMESPACE
BField::BField(ParamSet const &params)
: params{ params }
{
    this->fill(params.geomtr.B0); // fill with background B
}

void BField::update(EField const &efield, Real const dt) noexcept
{
    Real const cdtODx = dt * params.c / params.Dx;
    impl_update(*this, efield, cdtODx);
}

void BField::impl_update(BField &B, EField const &E, Real const cdtODx) noexcept
{
    auto const curl_E_times_cdt = [cdtODx](Vector const &E1, Vector const &E0) noexcept -> Vector {
        return {
            0,
            (-E1.z + E0.z) * cdtODx,
            (+E1.y - E0.y) * cdtODx,
        };
    };
    for (long i = 0; i < B.size(); ++i) {
        B[i] -= curl_E_times_cdt(E[i + 1], E[i + 0]);
    }
}

namespace {
template <class Object>
decltype(auto) write_attr(Object &obj, [[maybe_unused]] BField const &bfield)
{
    return obj;
}
} // namespace
auto operator<<(hdf5::Group &obj, [[maybe_unused]] BField const &bfield) -> decltype(obj)
{
    return write_attr(obj, bfield);
}
auto operator<<(hdf5::Dataset &obj, [[maybe_unused]] BField const &bfield) -> decltype(obj)
{
    return write_attr(obj, bfield);
}
PIC1D_END_NAMESPACE
