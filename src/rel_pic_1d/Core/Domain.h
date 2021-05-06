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

#ifndef Domain_h
#define Domain_h

#include "../Geometry.h"
#include "../ParamSet.h"
#include "./BField.h"
#include "./ColdSpecies.h"
#include "./Current.h"
#include "./EField.h"
#include "./PartSpecies.h"
#include "./Species.h"

#include <array>
#include <tuple>
#include <utility>

PIC1D_BEGIN_NAMESPACE
class Delegate;

class Domain {
public:
    ParamSet const                                          params;
    Geometry const                                          geomtr;
    Delegate *const                                         delegate;
    BField                                                  bfield;
    EField                                                  efield;
    Current                                                 current;
    std::array<PartSpecies, ParamSet::part_indices::size()> part_species;
    std::array<ColdSpecies, ParamSet::cold_indices::size()> cold_species;

private:
    BField  bfield_1; // temporary B at half time step
    Current J;
    bool    is_recurring_pass{ false };

public:
    ~Domain();
    explicit Domain(ParamSet const &params, Delegate *delegate);

    void advance_by(unsigned n_steps);

private:
    void cycle(Domain const &domain);

    template <class Species> Current const &collect_smooth(Current &J, Species const &sp) const;

    template <class... Ts, class Int, Int... Is>
    static auto make_part_species(ParamSet const &params, std::tuple<Ts...> const &descs,
                                  std::integer_sequence<Int, Is...>);
    template <class... Ts, class Int, Int... Is>
    static auto make_cold_species(ParamSet const &params, std::tuple<Ts...> const &descs,
                                  std::integer_sequence<Int, Is...>);
};
PIC1D_END_NAMESPACE

#endif /* Domain_h */
