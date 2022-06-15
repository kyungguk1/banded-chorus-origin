/*
 * Copyright (c) 2019-2022, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "EnergyRecorder.h"
#include <PIC/UTL/println.h>

#include <cmath>
#include <filesystem>
#include <functional>
#include <stdexcept>

PIC1D_BEGIN_NAMESPACE
auto EnergyRecorder::filepath(std::string_view const &wd) const
{
    constexpr std::string_view filename = "energy.csv";
    if (is_world_master())
        return std::filesystem::path{ wd } / filename;
    return std::filesystem::path{ "/dev/null" };
}

EnergyRecorder::EnergyRecorder(parallel::mpi::Comm _subdomain_comm, parallel::mpi::Comm const &world_comm, ParamSet const &params)
: Recorder{ Input::energy_recording_frequency, std::move(_subdomain_comm), world_comm }
{
    // open output stream
    //
    auto const path = filepath(params.working_directory);
    if (os.open(path, params.snapshot_load ? os.app : os.trunc); !os)
        throw std::invalid_argument{ std::string{ __PRETTY_FUNCTION__ } + " - open failed: " + path.string() };

    os.setf(os.scientific);
    os.precision(15);

    if (!params.snapshot_load) {
        // header lines
        //
        print(os, "step");   // integral step count
        print(os, ", time"); // simulation time
        //
        print(os, ", dB1^2/2, dB2^2/2, dB3^2/2"); // spatial average of fluctuating (without background) magnetic field energy density
        print(os, ", dE1^2/2, dE2^2/2, dE3^2/2"); // spatial average of fluctuating (without background) electric field energy density
        //
        for (unsigned i = 1; i <= ParamSet::part_indices::size(); ++i) {
            // spatial average of i'th species kinetic energy density in lab frame
            print(os, ", part_species(", i, ") KEt_lab", ", part_species(", i, ") KE1_lab", ", part_species(", i, ") KE2_lab", ", part_species(", i, ") KE3_lab");
            // spatial average of i'th species kinetic energy density in plasma frame
            print(os, ", part_species(", i, ") KEt_pla", ", part_species(", i, ") KE1_pla", ", part_species(", i, ") KE2_pla", ", part_species(", i, ") KE3_pla");
        }
        //
        for (unsigned i = 1; i <= ParamSet::cold_indices::size(); ++i) {
            // spatial average of i'th species kinetic energy density in lab frame
            print(os, ", cold_species(", i, ") KEt_lab", ", cold_species(", i, ") KE1_lab", ", cold_species(", i, ") KE2_lab", ", cold_species(", i, ") KE3_lab");
            // spatial average of i'th species kinetic energy density in plasma frame
            print(os, ", cold_species(", i, ") KEt_pla", ", cold_species(", i, ") KE1_pla", ", cold_species(", i, ") KE2_pla", ", cold_species(", i, ") KE3_pla");
        }
        //
        os << std::endl;
    }
}

void EnergyRecorder::record(const Domain &domain, const long step_count)
{
    if (!should_record_at(step_count))
        return;

    print(os, step_count, ", ", step_count * domain.params.dt);

    auto printer = [&os = this->os](MFAVector const &v) {
        print(os, ", ", v.x, ", ", v.y, ", ", v.z);
    };
    auto const &comm = subdomain_comm;

    printer(*comm.all_reduce(std::plus{}, dump(domain.bfield)));
    printer(*comm.all_reduce(std::plus{}, dump(domain.efield)));

    for (Species const &sp : domain.part_species) {
        auto const [KE_lab, KE_pla]
            = comm.all_reduce(std::plus{}, dump(sp)).unpack([](std::vector<FourMFATensor> vec) {
                  return std::make_pair(vec.at(0), vec.at(1));
              });
        print(os, ", ", KE_lab.tt); // total in lab frame
        printer(KE_lab.ts);         // kinetic in lab frame
        print(os, ", ", KE_pla.tt); // total in plasma frame
        printer(KE_pla.ts);         // kinetic in plasma frame
    }

    for (Species const &sp : domain.cold_species) {
        auto const [KE_lab, KE_pla]
            = comm.all_reduce(std::plus{}, dump(sp)).unpack([](std::vector<FourMFATensor> vec) {
                  return std::make_pair(vec.at(0), vec.at(1));
              });
        print(os, ", ", KE_lab.tt); // total in lab frame
        printer(KE_lab.ts);         // kinetic in lab frame
        print(os, ", ", KE_pla.tt); // total in plasma frame
        printer(KE_pla.ts);         // kinetic in plasma frame
    }

    os << std::endl;
}

auto EnergyRecorder::dump(BField const &bfield) -> MFAVector
{
    MFAVector  dB2O2{};
    auto const q1min = bfield.grid_subdomain_extent().min();
    for (long i = 0; i < bfield.size(); ++i) {
        auto const dB = bfield.geomtr.cart_to_mfa(bfield[i], CurviCoord{ i + q1min });
        dB2O2 += dB * dB;
    }
    dB2O2 /= 2 * Input::Nx;
    return dB2O2;
}
auto EnergyRecorder::dump(EField const &efield) -> MFAVector
{
    MFAVector  dE2O2{};
    auto const q1min = efield.grid_subdomain_extent().min();
    for (long i = 0; i < efield.size(); ++i) {
        auto const dE = efield.geomtr.cart_to_mfa(efield[i], CurviCoord{ i + q1min });
        dE2O2 += dE * dE;
    }
    dE2O2 /= 2 * Input::Nx;
    return dE2O2;
}
auto EnergyRecorder::dump(Species const &sp) -> std::vector<FourMFATensor>
{
    FourMFATensor KE_lab{};
    FourMFATensor KE_pla{};
    auto const    q1min = sp.grid_subdomain_extent().min();
    for (long i = 0; i < sp.moment<0>().size(); ++i) {
        auto const pos = CurviCoord{ i + q1min };
        auto const n   = Real{ sp.moment<0>()[i] };
        if (constexpr auto zero = 1e-15; n < zero)
            continue;

        auto const Mij_lab = sp.geomtr.cart_to_mfa(sp.moment<2>()[i], pos);
        {
            auto const E0 = n * sp.params.c2;
            KE_lab.tt += Mij_lab.tt - E0;
            KE_lab.ts += Mij_lab.ss.lo() / 2;
            KE_lab.ss += Mij_lab.ss;
        }

        auto const VOc     = sp.geomtr.cart_to_mfa(sp.moment<1>()[i], pos) / (n * sp.params.c);
        auto const beta    = std::sqrt(dot(VOc, VOc));
        auto const gamma   = 1 / std::sqrt((1 - beta) * (1 + beta));
        auto const Mij_pla = lorentz_boost<+1>(Mij_lab, VOc, gamma);
        {
            auto const E0 = n / gamma * sp.params.c2;
            KE_pla.tt += Mij_pla.tt - E0;
            KE_pla.ts += Mij_pla.ss.lo() / 2;
            KE_pla.ss += Mij_pla.ss;
        }
    }
    KE_lab *= sp.energy_density_conversion_factor() / Input::Nx;
    KE_pla *= sp.energy_density_conversion_factor() / Input::Nx;
    return { KE_lab, KE_pla };
}
PIC1D_END_NAMESPACE
