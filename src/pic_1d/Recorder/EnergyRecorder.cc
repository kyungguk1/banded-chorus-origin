/*
 * Copyright (c) 2019-2021, Kyungguk Min
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

#include "EnergyRecorder.h"

#include "../Utility/println.h"

#include <stdexcept>

// MARK:- P1D::EnergyRecorder
//
std::string P1D::EnergyRecorder::filepath(std::string const &wd) const
{
    constexpr char filename[] = "energy.csv";
    return is_master() ? wd + "/" + filename : null_dev;
}

P1D::EnergyRecorder::EnergyRecorder(parallel::mpi::Comm _comm, ParamSet const &params)
: Recorder{Input::energy_recording_frequency, std::move(_comm)}
{
    // open output stream
    //
    std::string const path = filepath(params.working_directory);
    if (os.open(path, params.snapshot_load ? os.app : os.trunc); !os)
        throw std::invalid_argument{std::string{__PRETTY_FUNCTION__} + " - open failed: " + path};

    os.setf(os.scientific);
    os.precision(15);

    if (!params.snapshot_load) {
        // header lines
        //
        print(os, "step");   // integral step count
        print(os, ", time"); // simulation time
        //
        print(os, ", dB1^2/2, dB2^2/2, dB3^2/2"); // spatial average of fluctuating (without
        // background) magnetic field energy density
        print(os, ", dE1^2/2, dE2^2/2, dE3^2/2"); // spatial average of fluctuating (without
        // background) electric field energy density
        //
        for (unsigned i = 1; i <= ParamSet::part_indices::size(); ++i) {
            // spatial average of i'th species kinetic energy density
            print(os, ", part_species(", i, ") mv1^2/2", ", part_species(", i, ") mv2^2/2",
                  ", part_species(", i, ") mv3^2/2");
            // spatial average of i'th species bulk flow energy density
            print(os, ", part_species(", i, ") mU1^2/2", ", part_species(", i, ") mU2^2/2",
                  ", part_species(", i, ") mU3^2/2");
        }
        //
        for (unsigned i = 1; i <= ParamSet::cold_indices::size(); ++i) {
            // spatial average of i'th species kinetic energy density
            print(os, ", cold_species(", i, ") mv1^2/2", ", cold_species(", i, ") mv2^2/2",
                  ", cold_species(", i, ") mv3^2/2");
            // spatial average of i'th species bulk flow energy density
            print(os, ", cold_species(", i, ") mU1^2/2", ", cold_species(", i, ") mU2^2/2",
                  ", cold_species(", i, ") mU3^2/2");
        }
        //
        os << std::endl;
    }
}

void P1D::EnergyRecorder::record(const Domain &domain, const long step_count)
{
    if (step_count % recording_frequency)
        return;

    print(os, step_count, ", ", step_count * domain.params.dt);

    auto printer = [&os = this->os](Vector const &v) {
        print(os, ", ", v.x, ", ", v.y, ", ", v.z);
    };
    using parallel::mpi::ReduceOp;

    printer(*comm.all_reduce<Vector>(ReduceOp::plus<Vector>(true), dump(domain.bfield)));
    printer(*comm.all_reduce<Vector>(ReduceOp::plus<Vector>(true), dump(domain.efield)));

    for (Species const &sp : domain.part_species) {
        Tensor const t = comm.all_reduce(ReduceOp::plus<Tensor>(true), dump(sp));
        printer(t.lo()); // kinetic
        printer(t.hi()); // bulk flow
    }

    for (Species const &sp : domain.cold_species) {
        Tensor const t = comm.all_reduce(ReduceOp::plus<Tensor>(true), dump(sp));
        printer(t.lo()); // kinetic
        printer(t.hi()); // bulk flow
    }

    os << std::endl;
}

auto P1D::EnergyRecorder::dump(BField const &bfield) noexcept -> Vector
{
    Vector dB2O2{};
    for (Vector const &B_ : bfield) {
        Vector const dB = bfield.geomtr.cart2fac(B_ - bfield.geomtr.B0);
        dB2O2 += dB * dB;
    }
    dB2O2 /= 2 * Input::Nx;
    return dB2O2;
}
auto P1D::EnergyRecorder::dump(EField const &efield) noexcept -> Vector
{
    Vector dE2O2{};
    for (Vector const &E_ : efield) {
        Vector const dE = efield.geomtr.cart2fac(E_);
        dE2O2 += dE * dE;
    }
    dE2O2 /= 2 * Input::Nx;
    return dE2O2;
}
auto P1D::EnergyRecorder::dump(Species const &sp) noexcept -> Tensor
{
    Tensor  KE{};
    Vector &mv2O2 = KE.lo(), &mU2O2 = KE.hi();
    for (long i = 0; i < sp.moment<0>().size(); ++i) {
        Real const     n{sp.moment<0>()[i]};
        Vector const   nV   = sp.geomtr.cart2fac(sp.moment<1>()[i]);
        Vector const   nvv  = sp.geomtr.cart2fac(sp.moment<2>()[i]);
        constexpr Real zero = 1e-15;
        mU2O2 += nV * nV / (n > zero ? n : 1);
        mv2O2 += nvv;
    }
    KE *= sp.energy_density_conversion_factor() / (2 * Input::Nx);
    return KE;
}
