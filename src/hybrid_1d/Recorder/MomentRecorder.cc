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

#include "MomentRecorder.h"

#include "../Utility/println.h"

#include <stdexcept>

std::string H1D::MomentRecorder::filepath(std::string const &wd, long const step_count) const
{
    constexpr char    prefix[] = "moment";
    std::string const filename = std::string{prefix} + "-" + std::to_string(step_count) + ".csv";
    return is_master() ? wd + "/" + filename : null_dev;
}

H1D::MomentRecorder::MomentRecorder(unsigned const rank, unsigned const size)
: Recorder{Input::moment_recording_frequency, rank, size}
{
    // configure output stream
    //
    os.setf(os.scientific);
    os.precision(15);
}

void H1D::MomentRecorder::record(const Domain &domain, const long step_count)
{
    if (step_count % recording_frequency)
        return;
    //
    std::string const path = filepath(domain.params.working_directory, step_count);
    if (os.open(path, os.trunc); !os) {
        throw std::invalid_argument{std::string{__FUNCTION__} + " - open failed: " + path};
    } else {
        // header lines
        //
        print(os, "step = ", step_count, "; ");
        print(os, "time = ", step_count * domain.params.dt, "; ");
        print(os, "Dx = ", domain.params.Dx, "; ");
        print(os, "Nx = ", domain.params.Nx, "; ");
        print(os, "Ns = ", domain.part_species.size() + domain.cold_species.size(), '\n');
        //
        for (unsigned i = 1; i <= domain.part_species.size(); ++i) {
            if (i - 1) {
                print(os, ", ");
            }
            print(os, "part_species(", i, ") <1>");
            print(os, ", part_species(", i, ") <v1>", ", part_species(", i, ") <v2>",
                  ", part_species(", i, ") <v3>");
            print(os, ", part_species(", i, ") <v1v1>", ", part_species(", i, ") <v2v2>",
                  ", part_species(", i, ") <v3v3>");
        }
        if (!domain.part_species.empty() && !domain.cold_species.empty()) {
            print(os, ", ");
        }
        for (unsigned i = 1; i <= domain.cold_species.size(); ++i) {
            if (i - 1) {
                print(os, ", ");
            }
            print(os, "cold_species(", i, ") <1>");
            print(os, ", cold_species(", i, ") <v1>", ", cold_species(", i, ") <v2>",
                  ", cold_species(", i, ") <v3>");
            print(os, ", cold_species(", i, ") <v1v1>", ", cold_species(", i, ") <v2v2>",
                  ", cold_species(", i, ") <v3v3>");
        }
        //
        print(os, '\n');

        // contents
        //
        auto printer = [&os = this->os](Vector const &v) -> std::ostream & {
            return print(os, v.x, ", ", v.y, ", ", v.z);
        };
        //
        auto tk = comm.send(
            std::make_pair(domain.part_species.begin(), domain.cold_species.begin()), master);
        if (is_master()) {
            using Payload = std::pair<PartSpecies const *, ColdSpecies const *>;
            comm.for_each<Payload>(
                all_ranks,
                [&os = this->os, Nx = domain.bfield.size(), Ns_part = domain.part_species.size(),
                 Ns_cold = domain.cold_species.size()](Payload payload, auto printer) {
                    auto [part_species, cold_species] = payload;
                    for (long i = 0; i < Nx; ++i) {
                        for (unsigned s = 0; s < Ns_part; ++s) {
                            if (s) {
                                print(os, ", ");
                            }
                            Species const &sp = part_species[s];
                            print(os, Real{sp.moment<0>()[i]}, ", ");
                            printer(sp.geomtr.cart2fac(sp.moment<1>()[i])) << ", ";
                            printer(sp.geomtr.cart2fac(sp.moment<2>()[i]));
                        }
                        if (Ns_part > 0 && Ns_cold > 0) {
                            print(os, ", ");
                        }
                        for (unsigned s = 0; s < Ns_cold; ++s) {
                            if (s) {
                                print(os, ", ");
                            }
                            Species const &sp = cold_species[s];
                            print(os, Real{sp.moment<0>()[i]}, ", ");
                            printer(sp.geomtr.cart2fac(sp.moment<1>()[i])) << ", ";
                            printer(sp.geomtr.cart2fac(sp.moment<2>()[i]));
                        }
                        //
                        print(os, '\n');
                    }
                },
                printer);
        }
        std::move(tk).wait();
    }
    os.close();
}
