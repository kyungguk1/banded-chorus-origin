/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "Domain.h"
#include "../Boundary/Delegate.h"

#include <cmath>

PIC1D_BEGIN_NAMESPACE
namespace {
template <class T, long N>
auto &operator+=(Grid<T, N, Pad> &lhs, Grid<T, N, Pad> const &rhs) noexcept
{
    auto rhs_first = rhs.dead_begin(), rhs_last = rhs.dead_end();
    auto lhs_first = lhs.dead_begin();
    while (rhs_first != rhs_last) {
        *lhs_first++ += *rhs_first++;
    }
    return lhs;
}
//
template <class T, long N> auto &operator*=(Grid<T, N, Pad> &lhs, T const rhs) noexcept
{
    auto first = lhs.dead_begin(), last = lhs.dead_end();
    while (first != last) {
        *first++ *= rhs;
    }
    return lhs;
}
} // namespace

// Domain impl
//
Domain::~Domain()
{
}
template <class... Ts, class Int, Int... Is>
auto Domain::make_part_species(ParamSet const &params, Geometry const &geo,
                               std::tuple<Ts...> const &descs, std::integer_sequence<Int, Is...>)
{
    static_assert((... && std::is_base_of_v<KineticPlasmaDesc, Ts>));
    static_assert(sizeof...(Ts) == sizeof...(Is));
    //
    auto const extent = params.Nx * Range{ 0, 1 };
    return std::array<PartSpecies, sizeof...(Ts)>{
        PartSpecies{ params, std::get<Is>(descs),
                     VDFVariant::make(std::get<Is>(descs), geo, extent, params.c) }...,
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
, geomtr{ params.O0, params.theta * M_PI / 180 }
, delegate{ delegate }
, bfield{ params }
, efield{ params }
, current{ params }
, part_species{}
, cold_species{}
, bfield_1{ params }
, J{ params }
{
    part_species = make_part_species(params, geomtr, params.part_descs, ParamSet::part_indices{});
    cold_species = make_cold_species(params, params.cold_descs, ParamSet::cold_indices{});
}

void Domain::advance_by(unsigned const n_steps)
{
    Domain &domain = *this;

    // pre-process
    //
    if (!is_recurring_pass) { // execute only once
        is_recurring_pass = true;
        delegate->once(domain);

        // fill in ghost cells
        //
        delegate->pass(domain, efield);
        delegate->pass(domain, bfield);
        for (ColdSpecies &sp : cold_species) {
            delegate->pass(domain, sp);
        }
    }

    // cycle
    //
    for (long i = 0; i < n_steps; ++i) {
        delegate->prologue(domain, i);
        cycle(domain);
        delegate->epilogue(domain, i);
    }

    // post-process; collect all moments
    //
    for (PartSpecies &sp : part_species) {
        sp.collect_all();
        delegate->gather(domain, sp);
    }
    for (ColdSpecies &sp : cold_species) {
        sp.collect_all();
    }
}
void Domain::cycle(Domain const &domain)
{
    BField     &bfield_0 = bfield;
    Real const &dt       = params.dt;
    //
    // 1. update B<0> from n-1/2 to n+1/2 using E at n
    //    B<1> = (B(n-1/2) + B(n+1/2))/2, so B<1> is at a full time step (n)
    //
    bfield_1 = bfield_0;
    {
        bfield_0.update(efield, dt);
        delegate->pass(domain, bfield_0);
    }
    (bfield_1 += bfield_0) *= Vector{ .5 };
    //
    // 2 & 3. update velocities and positions by dt and collect current density
    //
    current.reset();
    for (PartSpecies &sp : part_species) {
        sp.update_vel(bfield_1, efield, dt); // v(n-1/2) -> v(n+1/2)

        sp.update_pos(0.5 * dt, 0.5); // x(n) -> x(n+1/2)
        delegate->pass(domain, sp);

        sp.collect_part();
        current += collect_smooth(J, sp); // J(n+1/2)

        sp.update_pos(0.5 * dt, 0.5); // x(n+1/2) -> x(n+1)
        delegate->pass(domain, sp);
    }
    for (ColdSpecies &sp : cold_species) {
        sp.update_vel(bfield_1, efield, dt); // <v>(n-1/2) -> <v>(n+1/2)

        sp.collect_part();
        current += collect_smooth(J, sp); // J(n+1/2)
    }
    //
    // 5. update E from n to n+1 using B and J at n+1/2
    //
    efield.update(bfield_0, current, dt);
    delegate->pass(domain, efield);
}
template <class Species>
auto Domain::collect_smooth(Current &J, Species const &sp) const -> Current const &
{
    J.reset();
    //
    // collect & gather J
    //
    J += sp;
    delegate->gather(*this, J);
    //
    // optional smoothing
    //
    for (long i = 0; i < sp->number_of_source_smoothings; ++i) {
        delegate->pass(*this, J);
        J.smooth();
    }
    //
    delegate->pass(*this, J);
    return J;
}
PIC1D_END_NAMESPACE
