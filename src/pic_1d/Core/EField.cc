/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "EField.h"
#include "BField.h"
#include "Current.h"

PIC1D_BEGIN_NAMESPACE
EField::EField(ParamSet const &params)
: params{ params }
{
}

void EField::update(BField const &bfield, Current const &current, Real const dt) noexcept
{
    Real const cdtODx = dt * params.c / params.Dx;
    impl_update(*this, bfield, cdtODx, current, dt);
}

void EField::impl_update(EField &E, BField const &B, Real const cdtODx, Current const &J, Real const dt) noexcept
{
    auto const curl_B_times_cdt = [cdtODx](Vector const &B1, Vector const &B0) noexcept -> Vector {
        return {
            0,
            (-B1.z + B0.z) * cdtODx,
            (+B1.y - B0.y) * cdtODx,
        };
    };
    for (long i = 0; i < E.size(); ++i) {
        E[i] += curl_B_times_cdt(B[i - 0], B[i - 1]);
        E[i] -= J[i] * dt;
    }
}

namespace {
template <class Object>
decltype(auto) write_attr(Object &obj, [[maybe_unused]] EField const &efield)
{
    return obj;
}
} // namespace
auto operator<<(hdf5::Group &obj, [[maybe_unused]] EField const &efield) -> decltype(obj)
{
    return write_attr(obj, efield);
}
auto operator<<(hdf5::Dataset &obj, [[maybe_unused]] EField const &efield) -> decltype(obj)
{
    return write_attr(obj, efield);
}
PIC1D_END_NAMESPACE
