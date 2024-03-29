/*
 * Copyright (c) 2019-2022, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "Current.h"
#include "Species.h"

PIC1D_BEGIN_NAMESPACE
namespace {
template <class LIt, class RIt, class U>
void accumulate(LIt lhs_first, RIt rhs_first, RIt const rhs_last, U const &weight) noexcept
{
    while (rhs_first != rhs_last) {
        *lhs_first++ += *rhs_first++ * weight;
    }
}
} // namespace

Current::Current(ParamSet const &params)
: params{ params }
{
}

// current collector
//
auto Current::operator+=(Species const &sp) noexcept -> Current &
{
    accumulate(this->dead_begin(), sp.moment<1>().dead_begin(), sp.moment<1>().dead_end(), sp.current_density_conversion_factor());
    return *this;
}
PIC1D_END_NAMESPACE
