/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include "../ParamSet.h"
#include "BField.h"
#include "ColdSpecies.h"
#include "Current.h"
#include "EField.h"
#include "PartSpecies.h"
#include "Species.h"

#include <array>
#include <tuple>
#include <utility>

PIC1D_BEGIN_NAMESPACE
class Delegate;

class Domain {
public:
    ParamSet const                                          params;
    Delegate *const                                         delegate;
    BField                                                  bfield;
    EField                                                  efield;
    Current                                                 current;
    std::array<PartSpecies, ParamSet::part_indices::size()> part_species;
    std::array<ColdSpecies, ParamSet::cold_indices::size()> cold_species;

private:
    BField  bfield_1; // temporary B at half the time step
    Current J;
    bool    is_recurring_pass{ false };

public:
    ~Domain();
    Domain(ParamSet const &params, Delegate *delegate);

    void advance_by(unsigned n_steps);

private:
    void cycle(Domain const &domain);

    template <class Species>
    Current const &collect_smooth(Current &J, Species const &sp) const;

    template <class... Ts, class Int, Int... Is>
    static auto make_part_species(ParamSet const &params, std::tuple<Ts...> const &descs, std::integer_sequence<Int, Is...>);
    template <class... Ts, class Int, Int... Is>
    static auto make_cold_species(ParamSet const &params, std::tuple<Ts...> const &descs, std::integer_sequence<Int, Is...>);
};
PIC1D_END_NAMESPACE
