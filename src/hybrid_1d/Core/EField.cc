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
EField::EField(ParamSet const &params)
: params{ params }, geomtr{ params.geomtr }
{
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
auto EField::cart_to_contr(VectorGrid &B_contr, BField const &B_cart) const noexcept -> VectorGrid &
{
    constexpr auto ghost_offset = 1;
    static_assert(ghost_offset <= Pad);
    auto const q1min = params.full_grid_subdomain_extent.min();
    for (long i = -ghost_offset; i < BField::size() + ghost_offset; ++i) {
        B_contr[i] = geomtr.cart_to_contr(B_cart[i], CurviCoord{ i + q1min });
    }
    return B_contr;
}

void EField::update(BField const &bfield, Charge const &charge, Current const &current) noexcept
{
    impl_update_Pe(Pe, charge);
    mask(Pe, params.phase_retardation);

    impl_update_Je(Je, current, cart_to_covar(buffer, bfield));
    mask(Je, params.phase_retardation);

    impl_update_E(*this, Je, charge, cart_to_contr(buffer, bfield));
    mask(*this, params.amplitude_damping);
}
template <class T, long N>
void EField::mask(Grid<T, N, Pad> &E, MaskingFunction const &masking_function) const
{
    auto const left_offset  = params.half_grid_subdomain_extent.min() - params.half_grid_whole_domain_extent.min();
    auto const right_offset = params.half_grid_whole_domain_extent.max() - params.half_grid_subdomain_extent.max();
    for (long i = 0, first = 0, last = EField::size() - 1; i < EField::size(); ++i) {
        E[first++] *= masking_function(left_offset + i);
        E[last--] *= masking_function(right_offset + i);
    }
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
void EField::impl_update_Je(VectorGrid &Je_contr, Current const &Ji_cart, VectorGrid const &B_covar) const noexcept
{
    auto const curl_B_times_c
        = [cOsqrtg = params.c / geomtr.sqrt_g()](Vector const &B1, Vector const &B0) noexcept -> Vector {
        return {
            0,
            (-B1.z + B0.z) * cOsqrtg,
            (+B1.y - B0.y) * cOsqrtg,
        };
    };
    auto const q1min = params.half_grid_subdomain_extent.min();
    for (long i = 0; i < EField::size(); ++i) {
        Je_contr[i] = curl_B_times_c(B_covar[i + 1], B_covar[i + 0])
                    - geomtr.cart_to_contr(Ji_cart[i], CurviCoord{ i + q1min });
    }
}
void EField::impl_update_E(EField &E_cart, VectorGrid const &Je_contr, Charge const &rho, VectorGrid const &dB_contr) const noexcept
{
    auto const grad_Pe_times_c
        = [cO2 = params.c / 2](Scalar const &Pe_p1, Scalar const &Pe_m1) noexcept -> Vector {
        return {
            Real{ Pe_p1 - Pe_m1 } * cO2,
            0,
            0,
        };
    };
    auto const q1min = params.half_grid_subdomain_extent.min();
    for (long i = 0; i < EField::size(); ++i) {
        CurviCoord const pos{ i + q1min };

        // 1. Je x B term
        //
        auto const B_contr = geomtr.Bcontr(pos) + (dB_contr[i + 1] + dB_contr[i + 0]) * 0.5;
        auto const lorentz = cross(Je_contr[i], B_contr) * geomtr.sqrt_g();

        // 2. pressure gradient term
        //
        auto const pressure_force = grad_Pe_times_c(Pe[i - 1], Pe[i + 1]); // reverse order is to account for pressure force = -grad P

        // 3. electric field
        //
        auto const E_covar = (lorentz + pressure_force) / Real{ rho[i] };

        E_cart[i] = geomtr.covar_to_cart(E_covar, pos);
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
