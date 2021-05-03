/*
 * Copyright (c) 2019, Kyungguk Min
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * - Redistributions of source code must retain the above copyright notice, this
 *  list of conditions and the following disclaimer.
 *
 * - Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

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
        PartSpecies{ params, std::get<Is>(descs), VDF::make(std::get<Is>(descs)) }...,
    };
}
template <class... Ts, class Int, Int... Is>
auto H1D::Domain::make_cold_species(ParamSet const &params, std::tuple<Ts...> const &descs,
                                    std::integer_sequence<Int, Is...>)
{
    static_assert((... && std::is_base_of_v<ColdPlasmaDesc, Ts>));
    static_assert(sizeof...(Ts) == sizeof...(Is));
    //
    return std::array<ColdSpecies, sizeof...(Ts)>{ ColdSpecies{ params, std::get<Is>(descs) }... };
}
H1D::Domain::Domain(ParamSet const &params, Delegate *delegate)
: params{ params }
, geomtr{ params }
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
