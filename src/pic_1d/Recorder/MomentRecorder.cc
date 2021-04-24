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

#include "MomentRecorder.h"

#include "../Utility/println.h"

#include <stdexcept>

// MARK:- P1D::MomentRecorder
//
std::string P1D::MomentRecorder::filepath(std::string const &wd, long const step_count) const
{
    constexpr char    prefix[] = "moment";
    std::string const filename = std::string{prefix} + "-" + std::to_string(step_count) + ".csv";
    return is_master() ? wd + "/" + filename : null_dev;
}

P1D::MomentRecorder::MomentRecorder(parallel::mpi::Comm _comm)
: Recorder{Input::moment_recording_frequency, std::move(_comm)}
{
    // configure output stream
    //
    os.setf(os.scientific);
    os.precision(15);
}

void P1D::MomentRecorder::record(const Domain &domain, const long step_count)
{
    if (step_count % recording_frequency)
        return;

    std::string const path = filepath(domain.params.working_directory, step_count);
    if (os.open(path, os.trunc); !os)
        throw std::invalid_argument{std::string{__PRETTY_FUNCTION__} + " - open failed: " + path};

    // header lines
    //
    print(os, "step = ", step_count, "; ");
    print(os, "time = ", step_count * domain.params.dt, "; ");
    print(os, "Dx = ", domain.params.Dx, "; ");
    print(os, "Nx = ", domain.params.Nx, "; ");
    print(os, "Ns = ", domain.part_species.size() + domain.cold_species.size(), '\n');

    for (unsigned long i = 1; i <= domain.part_species.size(); ++i) {
        if (i - 1)
            print(os, ", ");

        print(os, "part_species(", i, ") <1>");
        print(os, ", part_species(", i, ") <v1>", ", part_species(", i, ") <v2>", ", part_species(",
              i, ") <v3>");
        print(os, ", part_species(", i, ") <v1v1>", ", part_species(", i, ") <v2v2>",
              ", part_species(", i, ") <v3v3>");
    }

    if (!domain.part_species.empty() && !domain.cold_species.empty())
        print(os, ", ");

    for (unsigned long i = 1; i <= domain.cold_species.size(); ++i) {
        if (i - 1)
            print(os, ", ");

        print(os, "cold_species(", i, ") <1>");
        print(os, ", cold_species(", i, ") <v1>", ", cold_species(", i, ") <v2>", ", cold_species(",
              i, ") <v3>");
        print(os, ", cold_species(", i, ") <v1v1>", ", cold_species(", i, ") <v2v2>",
              ", cold_species(", i, ") <v3v3>");
    }

    print(os, '\n');

    // contents
    //
    auto printer = [&os = this->os](Vector const &v) -> std::ostream & {
        return print(os, v.x, ", ", v.y, ", ", v.z);
    };
    auto check_Nx = [pretty = std::string{__PRETTY_FUNCTION__},
                     Nx     = domain.params.Nx](auto const &v, char const *label) {
        if (v.size() != Nx)
            throw std::runtime_error{pretty + " - incorrect array size of " + label};
    };

    if (is_master()) {
        std::vector<std::vector<Scalar>> part_mom0, cold_mom0;
        std::vector<std::vector<Vector>> part_mom1, cold_mom1;
        std::vector<std::vector<Tensor>> part_mom2, cold_mom2;

        for (PartSpecies const &sp : domain.part_species) {
            check_Nx(part_mom0.emplace_back(
                         comm.gather<0>({sp.moment<0>().begin(), sp.moment<0>().end()}, master)),
                     "part_mom0");

            check_Nx(part_mom1.emplace_back(
                         comm.gather<1>({sp.moment<1>().begin(), sp.moment<1>().end()}, master)),
                     "part_mom1");

            check_Nx(part_mom2.emplace_back(
                         comm.gather<2>({sp.moment<2>().begin(), sp.moment<2>().end()}, master)),
                     "part_mom2");
        }

        for (ColdSpecies const &sp : domain.cold_species) {
            check_Nx(cold_mom0.emplace_back(
                         comm.gather<0>({sp.moment<0>().begin(), sp.moment<0>().end()}, master)),
                     "cold_mom0");

            check_Nx(cold_mom1.emplace_back(
                         comm.gather<1>({sp.moment<1>().begin(), sp.moment<1>().end()}, master)),
                     "cold_mom1");

            check_Nx(cold_mom2.emplace_back(
                         comm.gather<2>({sp.moment<2>().begin(), sp.moment<2>().end()}, master)),
                     "cold_mom2");
        }

        for (unsigned long i = 0; i < domain.params.Nx; ++i) {
            for (unsigned long s = 0; s < domain.part_species.size(); ++s) {
                if (s)
                    print(os, ", ");

                print(os, Real{part_mom0.at(s).at(i)}, ", ");
                printer(domain.geomtr.cart2fac(part_mom1.at(s).at(i))) << ", ";
                printer(domain.geomtr.cart2fac(part_mom2.at(s).at(i)));
            }

            if (!domain.part_species.empty() && !domain.cold_species.empty())
                print(os, ", ");

            for (unsigned long s = 0; s < domain.cold_species.size(); ++s) {
                if (s)
                    print(os, ", ");

                print(os, Real{cold_mom0.at(s).at(i)}, ", ");
                printer(domain.geomtr.cart2fac(cold_mom1.at(s).at(i))) << ", ";
                printer(domain.geomtr.cart2fac(cold_mom2.at(s).at(i)));
            }

            print(os, '\n');
        }
    } else {
        for (PartSpecies const &sp : domain.part_species) {
            comm.gather<0>(sp.moment<0>().begin(), sp.moment<0>().end(), nullptr, master);
            comm.gather<1>(sp.moment<1>().begin(), sp.moment<1>().end(), nullptr, master);
            comm.gather<2>(sp.moment<2>().begin(), sp.moment<2>().end(), nullptr, master);
        }
        for (ColdSpecies const &sp : domain.cold_species) {
            comm.gather<0>(sp.moment<0>().begin(), sp.moment<0>().end(), nullptr, master);
            comm.gather<1>(sp.moment<1>().begin(), sp.moment<1>().end(), nullptr, master);
            comm.gather<2>(sp.moment<2>().begin(), sp.moment<2>().end(), nullptr, master);
        }
    }

    os.close();
}
