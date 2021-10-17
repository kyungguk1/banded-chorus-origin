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
: params{ params }, geomtr{ params.geomtr }
{
}

void EField::update(BField const &bfield, Current const &current, Real const dt) noexcept
{
    Real const cdtOsqrtg = dt * params.c / geomtr.sqrt_g();
    impl_update(*this, cart_to_covar(buffer, bfield), cdtOsqrtg, current, dt);
}

void EField::impl_update(EField &E_cart, VectorGrid const &B_covar, Real const cdtOsqrtg, Current const &J_cart, Real const dt) const noexcept
{
    auto const curl_B_times_cdt = [cdtOsqrtg](Vector const &B1, Vector const &B0) noexcept -> Vector {
        return {
            0,
            (-B1.z + B0.z) * cdtOsqrtg,
            (+B1.y - B0.y) * cdtOsqrtg,
        };
    };
    auto const q1min = params.half_grid_subdomain_extent.min();
    for (long i = 0; i < EField::size(); ++i) {
        auto const E_contr = curl_B_times_cdt(B_covar[i + 1], B_covar[i + 0]);
        E_cart[i] += geomtr.contr_to_cart(E_contr, CurviCoord{ i + q1min });
        E_cart[i] -= J_cart[i] * dt;
    }
}
auto EField::cart_to_covar(VectorGrid &B_covar, BField const &B_cart) const noexcept -> VectorGrid &
{
    constexpr auto ghost_offset = 1;
    static_assert(ghost_offset <= Pad);
    auto const q1min = params.full_grid_subdomain_extent.min();
    for (long i = -ghost_offset; i < BField::size() + ghost_offset; ++i) {
        B_covar[i] = geomtr.cart_to_covar(B_cart[i], CurviCoord{ i + q1min });
    }
    return B_covar;
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
