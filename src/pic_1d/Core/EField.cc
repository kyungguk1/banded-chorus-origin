/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "EField.h"
#include "BField.h"
#include "Current.h"

#include <algorithm>
#include <functional>

PIC1D_BEGIN_NAMESPACE
EField::EField(ParamSet const &params)
: params{ params }, geomtr{ params.geomtr }
{
}

void EField::update(BField const &bfield, Current const &current, Real const dt) noexcept
{
    Real const cdtOsqrtg = dt * params.c / geomtr.sqrt_g();

    E_prev.operator=(*this);
    this->fill(Vector{});

    // Delta-E followed by phase retardation
    impl_update(*this, cart_to_covar(Bcovar, bfield), cdtOsqrtg, current, dt);
    mask(*this, params.phase_retardation);

    // Next-E followed by amplitude damping
    std::transform(this->begin(), this->end(), E_prev.begin(), this->begin(), std::plus{});
    mask(*this, params.amplitude_damping);
}
void EField::mask(EField &E, MaskingFunction const &masking_function) const
{
    auto const left_offset  = params.half_grid_subdomain_extent.min() - params.half_grid_whole_domain_extent.min();
    auto const right_offset = params.half_grid_whole_domain_extent.max() - params.half_grid_subdomain_extent.max();
    for (long i = 0, first = 0, last = EField::size() - 1; i < EField::size(); ++i) {
        E[first++] *= masking_function(left_offset + i);
        E[last--] *= masking_function(right_offset + i);
    }
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
