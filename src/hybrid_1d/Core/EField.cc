/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "EField.h"
#include "BField.h"
#include "Charge.h"
#include "Current.h"

#include <cmath>

HYBRID1D_BEGIN_NAMESPACE
EField::EField(ParamSet const &params) : params{ params }
{
}

void EField::update(BField const &bfield, Charge const &charge, Current const &current) noexcept
{
    impl_update_Pe(Pe, charge);
    impl_update_Je(Je, current, bfield);
    impl_update_E(*this, bfield, charge);
}

void EField::impl_update_Pe(ScalarGrid &Pe, Charge const &rho) const noexcept
{
    eFluidDesc const ef = params.efluid_desc;
    //
    Real const O0        = params.O0;
    Real const O02beO2   = (O0 * O0) * ef.beta * 0.5;
    Real const mOeOO0oe2 = -ef.Oc / (O0 * (ef.op * ef.op));
    for (long i = -Pad; i < Pe.size() + Pad; ++i) {
        Pe[i] = std::pow(mOeOO0oe2 * Real{ rho[i] }, ef.gamma) * O02beO2;
    }
}
void EField::impl_update_Je(VectorGrid &Je, Current const &Ji, BField const &B) const noexcept
{
    Real const cODx = params.c / params.Dx;
    for (long i = 0; i < B.size(); ++i) {
        // J total
        //
        Je[i].x = 0;
        Je[i].y = (-B[i + 1].z + B[i + 0].z) * cODx;
        Je[i].z = (+B[i + 1].y - B[i + 0].y) * cODx;

        // Je = J - Ji
        //
        Je[i] -= Ji[i];
    }
}
void EField::impl_update_E(EField &E, BField const &B, Charge const &rho) const noexcept
{
    Real const cODx = params.c / params.Dx;
    for (long i = 0; i < E.size(); ++i) {
        Vector &Ei = E[i];

        // 1. Je x B term
        //
        Ei = cross(Je[i], (B[i + 1] + B[i + 0]) * 0.5);

        // 2. pressure gradient term
        //
        Ei.x -= 0.5 * Real{ Pe[i + 1] - Pe[i - 1] } * cODx;
        Ei.y -= 0;
        Ei.z -= 0;

        // 3. divide by charge density
        //
        Ei /= Real{ rho[i] };
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
HYBRID1D_END_NAMESPACE
