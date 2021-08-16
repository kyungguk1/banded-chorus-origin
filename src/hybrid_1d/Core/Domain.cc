/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "Domain.h"

HYBRID1D_BEGIN_NAMESPACE
Domain::~Domain()
{
}
template <class... Ts, class Int, Int... Is>
auto Domain::make_part_species(ParamSet const &params, std::tuple<Ts...> const &descs,
                               std::integer_sequence<Int, Is...>)
{
    static_assert((... && std::is_base_of_v<KineticPlasmaDesc, Ts>));
    static_assert(sizeof...(Ts) == sizeof...(Is));
    //
    auto const extent = params.Nx * Range{ 0, 1 };
    return std::array<PartSpecies, sizeof...(Ts)>{
        PartSpecies{ params, std::get<Is>(descs),
                     VDFVariant::make(std::get<Is>(descs), params.geomtr, extent, params.c) }...,
    };
}
template <class... Ts, class Int, Int... Is>
auto Domain::make_cold_species(ParamSet const &params, std::tuple<Ts...> const &descs,
                               std::integer_sequence<Int, Is...>)
{
    static_assert((... && std::is_base_of_v<ColdPlasmaDesc, Ts>));
    static_assert(sizeof...(Ts) == sizeof...(Is));
    //
    return std::array<ColdSpecies, sizeof...(Ts)>{ ColdSpecies{ params, std::get<Is>(descs) }... };
}
Domain::Domain(ParamSet const &params, Delegate *delegate)
: params{ params }
, delegate{ delegate }
, bfield{ params }
, efield{ params }
, charge{ params }
, current{ params }
, part_species{ make_part_species(params, params.part_descs, ParamSet::part_indices{}) }
, cold_species{ make_cold_species(params, params.cold_descs, ParamSet::cold_indices{}) }
, rho{ params }
, J{ params }
{
}
HYBRID1D_END_NAMESPACE
