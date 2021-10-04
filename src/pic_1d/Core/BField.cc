/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "BField.h"
#include "EField.h"

PIC1D_BEGIN_NAMESPACE
BField::BField(ParamSet const &params)
: params{ params }, geomtr{ params.geomtr }
{
}

auto BField::operator=(BField const &o) noexcept -> BField &
{ // this is only for the return type casting
    this->Grid::operator=(o);
    return *this;
}

void BField::update(EField const &efield, Real const dt) noexcept
{
    Real const cdtOsqrtg = dt * params.c / geomtr.sqrt_g();
    impl_update(*this, cart_to_covar(buffer, efield), cdtOsqrtg);
}

void BField::impl_update(BField &B_cart, VectorGrid const &E_covar, Real const cdtOsqrtg) const noexcept
{
    auto const curl_E_times_cdt = [cdtOsqrtg](Vector const &E1, Vector const &E0) noexcept -> Vector {
        return {
            0,
            (-E1.z + E0.z) * cdtOsqrtg,
            (+E1.y - E0.y) * cdtOsqrtg,
        };
    };
    auto const q1min = params.half_grid_subdomain_extent.min();
    for (long i = 0; i < BField::size(); ++i) {
        auto const B_contr = curl_E_times_cdt(E_covar[i + 1], E_covar[i + 0]);
        B_cart[i] -= geomtr.contr_to_cart(B_contr, CurviCoord{ i + q1min });
    }
}
auto BField::cart_to_covar(VectorGrid &E_covar, EField const &E_cart) const noexcept -> VectorGrid &
{
    constexpr auto ghost_offset = 1;
    static_assert(ghost_offset <= Pad);
    auto const q1min = params.full_grid_subdomain_extent.min();
    for (long i = -ghost_offset; i < EField::size() + ghost_offset; ++i) {
        E_covar[i] = geomtr.cart_to_covar(E_cart[i], CurviCoord{ i + q1min });
    }
    return E_covar;
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
