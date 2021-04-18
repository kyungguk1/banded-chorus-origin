//
// Copyright (c) 2019, Kyungguk Min
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// - Redistributions of source code must retain the above copyright notice, this
//  list of conditions and the following disclaimer.
//
// - Redistributions in binary form must reproduce the above copyright notice,
//  this list of conditions and the following disclaimer in the documentation
//  and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include "Domain_PC.h"

#include "./Domain.hh"

H1D::Domain_PC::Domain_PC(ParamSet const &params, Delegate *delegate)
: Domain{params, delegate}
, bfield_1{params}
, efield_1{params}
, part_predict{params}
, cold_predict{params}
{
}

void H1D::Domain_PC::advance_by(unsigned const n_steps)
{
    Domain &domain = *this;

    // pre-process
    //
    if (!is_recurring_pass) { // execute only once
        is_recurring_pass = true;
        delegate->once(domain);
        //
        // fill in ghost cells
        //
        delegate->pass(domain, efield);
        delegate->pass(domain, bfield);
        //
        // deposit charge and current densities
        //
        charge.reset();
        current.reset();
        for (PartSpecies &sp : part_species) {
            sp.collect_part();
            charge += collect_smooth(rho, sp);
            current += collect_smooth(J, sp);
        }
        for (ColdSpecies &sp : cold_species) {
            delegate->pass(domain, sp);
            sp.collect_part();
            charge += collect_smooth(rho, sp);
            current += collect_smooth(J, sp);
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
void H1D::Domain_PC::cycle(Domain const &domain)
{
    predictor_step(domain);
    corrector_step(domain);
}
void H1D::Domain_PC::predictor_step(Domain const &domain)
{
    BField &    bfield_0 = this->bfield;
    EField &    efield_0 = this->efield;
    Real const &dt       = params.dt;
    //
    // 1. Faraday's law; predict 1
    //
    bfield_1 = bfield_0;
    bfield_1.update(efield_0, dt), delegate->pass(domain, bfield_1);
    //
    // 2. Ohm's law; predict 1
    //
    efield_1.update(bfield_1, charge, current), delegate->pass(domain, efield_1);
    //
    // 3. Average fields
    //
    (bfield_1 += bfield_0) *= Vector{.5};
    (efield_1 += efield_0) *= Vector{.5};
    //
    // 4 & 5. Particle push and deposit charge and current densities; predict
    //
    charge.reset();
    current.reset();
    for (PartSpecies const &sp : part_species) {
        auto &predictor = part_predict = sp;

        predictor.update_pos(0.5 * dt, 0.5);
        predictor.update_vel(bfield_1, efield_1, dt);
        predictor.update_pos(0.5 * dt, 0.5);
        delegate->pass(domain, predictor);

        predictor.collect_part();
        charge += collect_smooth(rho, predictor);
        current += collect_smooth(J, predictor);
    }
    for (ColdSpecies const &sp : cold_species) {
        auto &predictor = cold_predict = sp;

        predictor.update_den(0.5 * dt);
        delegate->pass(domain, predictor);
        predictor.update_vel(bfield_1, efield_1, dt);
        delegate->pass(domain, predictor);
        predictor.update_den(0.5 * dt);
        delegate->pass(domain, predictor);

        predictor.collect_part();
        charge += collect_smooth(rho, predictor);
        current += collect_smooth(J, predictor);
    }
}
void H1D::Domain_PC::corrector_step(Domain const &domain)
{
    BField &    bfield_0 = this->bfield;
    EField &    efield_0 = this->efield;
    Real const &dt       = params.dt;
    //
    // 6. Faraday's law; predict 2
    //
    bfield_1 = bfield_0;
    bfield_1.update(efield_1, dt), delegate->pass(domain, bfield_1);
    //
    // 7. Ohm's law; predict 2
    //
    efield_1.update(bfield_1, charge, current), delegate->pass(domain, efield_1);
    //
    // 8. Average fields
    //
    (bfield_1 += bfield_0) *= Vector{.5};
    (efield_1 += efield_0) *= Vector{.5};
    //
    // 9 & 10. Particle push and deposit charge and current densities; correct
    //
    charge.reset();
    current.reset();
    for (PartSpecies &sp : part_species) {
        sp.update_pos(0.5 * dt, 0.5);
        sp.update_vel(bfield_1, efield_1, dt);
        sp.update_pos(0.5 * dt, 0.5);
        delegate->pass(domain, sp);

        sp.collect_part();
        charge += collect_smooth(rho, sp);
        current += collect_smooth(J, sp);
    }
    for (ColdSpecies &sp : cold_species) {
        sp.update_den(0.5 * dt);
        delegate->pass(domain, sp);
        sp.update_vel(bfield_1, efield_1, dt);
        delegate->pass(domain, sp);
        sp.update_den(0.5 * dt);
        delegate->pass(domain, sp);

        sp.collect_part();
        charge += collect_smooth(rho, sp);
        current += collect_smooth(J, sp);
    }
    //
    // 11. Faraday's law; correct
    //
    bfield_0.update(efield_1, dt), delegate->pass(domain, bfield_0);
    //
    // 12. Ohm's law; correct
    //
    efield_0.update(bfield_0, charge, current), delegate->pass(domain, efield_0);
}
