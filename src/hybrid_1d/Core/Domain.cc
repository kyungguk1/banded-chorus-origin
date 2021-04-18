//
//  Domain.cc
//  hybrid_1d
//
//  Created by KYUNGGUK MIN on 1/18/19.
//  Copyright © 2019 Kyungguk Min & Kaijun Liu. All rights reserved.
//

#include "Domain.h"

#include <cmath>

H1D::Domain::~Domain()
{
}
template <class... Ts, class Int, Int... Is>
auto H1D::Domain::make_part_species(ParamSet const &params, std::tuple<Ts...> const &descs,
                                    std::integer_sequence<Int, Is...>)
{
    static_assert((... && std::is_base_of_v<KineticPlasmaDesc, Ts>));
    static_assert(sizeof...(Ts) == sizeof...(Is));
    //
    return std::array<PartSpecies, sizeof...(Ts)>{
        PartSpecies{params, std::get<Is>(descs), VDF::make(std::get<Is>(descs))}...};
}
template <class... Ts, class Int, Int... Is>
auto H1D::Domain::make_cold_species(ParamSet const &params, std::tuple<Ts...> const &descs,
                                    std::integer_sequence<Int, Is...>)
{
    static_assert((... && std::is_base_of_v<ColdPlasmaDesc, Ts>));
    static_assert(sizeof...(Ts) == sizeof...(Is));
    //
    return std::array<ColdSpecies, sizeof...(Ts)>{ColdSpecies{params, std::get<Is>(descs)}...};
}
H1D::Domain::Domain(ParamSet const &params, Delegate *delegate)
: params{params}
, geomtr{params}
, delegate{delegate}
, bfield{params}
, efield{params}
, charge{params}
, current{params}
, part_species{make_part_species(params, params.part_descs, ParamSet::part_indices{})}
, cold_species{make_cold_species(params, params.cold_descs, ParamSet::cold_indices{})}
, rho{params}
, J{params}
{
}
