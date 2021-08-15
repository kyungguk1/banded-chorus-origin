/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "EnergyRecorder.h"
#include <PIC/println.h>

#include <stdexcept>

PIC1D_BEGIN_NAMESPACE
inline std::string EnergyRecorder::filepath(std::string const &wd) const
{
    constexpr char filename[] = "energy.csv";
    return is_master() ? wd + "/" + filename : null_dev;
}

EnergyRecorder::EnergyRecorder(parallel::mpi::Comm _comm, ParamSet const &params)
: Recorder{ Input::energy_recording_frequency, std::move(_comm) }
{
    // open output stream
    //
    std::string const path = filepath(params.working_directory);
    if (os.open(path, params.snapshot_load ? os.app : os.trunc); !os)
        throw std::invalid_argument{ std::string{ __PRETTY_FUNCTION__ }
                                     + " - open failed: " + path };

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

void EnergyRecorder::record(const Domain &domain, const long step_count)
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

auto EnergyRecorder::dump(BField const &bfield) noexcept -> Vector
{
    Vector dB2O2{};
    for (Vector const &B_ : bfield) {
        Vector const dB = bfield.params.geomtr.cart2fac(B_ - bfield.params.geomtr.B0);
        dB2O2 += dB * dB;
    }
    dB2O2 /= 2 * Input::Nx;
    return dB2O2;
}
auto EnergyRecorder::dump(EField const &efield) noexcept -> Vector
{
    Vector dE2O2{};
    for (Vector const &E_ : efield) {
        Vector const dE = efield.params.geomtr.cart2fac(E_);
        dE2O2 += dE * dE;
    }
    dE2O2 /= 2 * Input::Nx;
    return dE2O2;
}
auto EnergyRecorder::dump(Species const &sp) noexcept -> Tensor
{
    Tensor  KE{};
    Vector &mv2O2 = KE.lo(), &mU2O2 = KE.hi();
    for (long i = 0; i < sp.moment<0>().size(); ++i) {
        auto const n = Real{ sp.moment<0>()[i] };
        if (constexpr auto zero = 1e-15; n < zero)
            continue;

        Vector const nV  = sp.params.geomtr.cart2fac(sp.moment<1>()[i]);
        Vector const nvv = sp.params.geomtr.cart2fac(sp.moment<2>()[i]).lo();
        mU2O2 += nV * nV / n;
        mv2O2 += nvv;
    }
    KE *= sp.energy_density_conversion_factor() / (2 * Input::Nx);
    return KE;
}
PIC1D_END_NAMESPACE
