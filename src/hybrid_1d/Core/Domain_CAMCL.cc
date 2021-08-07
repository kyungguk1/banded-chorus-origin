/*
 * Copyright (c) 2019, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "Domain_CAMCL.h"

#include "./Domain.hh"

H1D::Domain_CAMCL::Domain_CAMCL(ParamSet const &params, Delegate *delegate)
: Domain{ params, delegate }
, bfield_1{ params }
, current_1{ params }
, charge_1{ params }
, lambda{ params }
, gamma{ params }
{
}

void H1D::Domain_CAMCL::advance_by(unsigned const n_steps)
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
void H1D::Domain_CAMCL::cycle(Domain const &domain)
{
    Current &   current_0 = this->current;
    Charge &    charge_0  = this->charge;
    Real const &dt        = params.dt;
    //
    // 1 & 2. update velocities and positions by full step and collect charge and current densities
    //
    current_0.reset();
    current_1.reset();
    charge_0.reset();
    charge_1.reset();
    for (PartSpecies &sp : part_species) {
        sp.update_vel(bfield, efield, dt); // v^N-1/2 -> v^N+1/2

        sp.collect_part();
        current_0 += collect_smooth(J, sp);  // J^-
        charge_0 += collect_smooth(rho, sp); // rho^N

        sp.update_pos(dt, 1); // x^N -> x^N+1
        delegate->pass(domain, sp);

        sp.collect_part();
        current_1 += collect_smooth(J, sp);  // J^+
        charge_1 += collect_smooth(rho, sp); // rho^N+1
    }
    for (ColdSpecies &sp : cold_species) {
        sp.update_vel(bfield, efield, dt); // v^N-1/2 -> v^N+1/2
        delegate->pass(domain, sp);

        sp.collect_part();
        current_0 += collect_smooth(J, sp);  // J^-
        charge_0 += collect_smooth(rho, sp); // rho^N

        sp.update_den(dt); // x^N -> x^N+1
        delegate->pass(domain, sp);

        sp.collect_part();
        current_1 += collect_smooth(J, sp);  // J^+
        charge_1 += collect_smooth(rho, sp); // rho^N+1
    }
    //
    // 3. average charge and current densities
    //
    (charge_0 += charge_1) *= Scalar{ .5 };   // rho^N+1/2
    (current_0 += current_1) *= Vector{ .5 }; // J^N+1/2
    //
    // 4. subcycle magnetic field by full step
    //
    subcycle(domain, charge_0, current_0, dt);
    //
    // 5. calculate electric field* and advance current density
    //
    efield.update(bfield, charge_1, current_0), delegate->pass(domain, efield);
    for (PartSpecies const &sp : part_species) {
        lambda.reset(), delegate->gather(domain, lambda += sp);
        gamma.reset(), delegate->gather(domain, gamma += sp);
        J.reset(), J.advance(lambda, gamma, bfield, efield, dt / 2.0);
        for (long i = 0; i < sp->number_of_source_smoothings; ++i) {
            delegate->pass(domain, J), J.smooth();
        }
        delegate->pass(domain, J), current_1 += J;
    }
    for (ColdSpecies const &sp : cold_species) {
        lambda.reset(), delegate->gather(domain, lambda += sp);
        gamma.reset(), delegate->gather(domain, gamma += sp);
        J.reset(), J.advance(lambda, gamma, bfield, efield, dt / 2.0);
        for (long i = 0; i < sp->number_of_source_smoothings; ++i) {
            delegate->pass(domain, J), J.smooth();
        }
        delegate->pass(domain, J), current_1 += J;
    }
    //
    // 6. calculate electric field
    //
    efield.update(bfield, charge_1, current_1), delegate->pass(domain, efield);
}
void H1D::Domain_CAMCL::subcycle(Domain const &domain, Charge const &charge, Current const &current,
                                 Real const _dt)
{
    BField &       bfield_0 = this->bfield;
    constexpr long m        = Input::n_subcycles;
    static_assert(m >= 2, "invalid n_subcycles");
    Real const dt = _dt / m, dt_x_2 = dt * 2.0;
    //
    // prologue
    //
    bfield_1 = bfield_0;
    efield.update(bfield_0, charge, current), delegate->pass(domain, efield);
    bfield_1.update(efield, dt), delegate->pass(domain, bfield_1);
    //
    // loop
    //
    for (long i = 1; i < m; ++i) {
        efield.update(bfield_1, charge, current), delegate->pass(domain, efield);
        bfield_0.update(efield, dt_x_2), delegate->pass(domain, bfield_0);
        bfield_0.swap(bfield_1);
    }
    //
    // epilogue
    //
    efield.update(bfield_1, charge, current), delegate->pass(domain, efield);
    bfield_0.update(efield, dt), delegate->pass(domain, bfield_0);
    //
    // average
    //
    (bfield_0 += bfield_1) *= Vector{ .5 };
}
