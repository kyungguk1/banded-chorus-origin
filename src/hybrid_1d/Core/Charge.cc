/*
 * Copyright (c) 2019, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "Charge.h"

#include "./Species.h"

// helper
//
namespace {
template <class LIt, class RIt, class U>
void accumulate(LIt lhs_first, RIt rhs_first, RIt const rhs_last, U const &weight) noexcept
{
    while (rhs_first != rhs_last) {
        *lhs_first++ += *rhs_first++ * weight;
    }
}
} // namespace

H1D::Charge::Charge(ParamSet const &params) : GridQ{}, tmp{}, params{ params }, geomtr{ params }
{
}

// density collector
//
H1D::Charge &H1D::Charge::operator+=(Species const &sp) noexcept
{
    ::accumulate(this->dead_begin(), sp.moment<0>().dead_begin(), sp.moment<0>().dead_end(),
                 sp.charge_density_conversion_factor());
    return *this;
}

H1D::Lambda &H1D::Lambda::operator+=(Species const &sp) noexcept
{
    ::accumulate(this->dead_begin(), sp.moment<0>().dead_begin(), sp.moment<0>().dead_end(),
                 sp.charge_density_conversion_factor() * sp->Oc / params.O0);
    return *this;
}
