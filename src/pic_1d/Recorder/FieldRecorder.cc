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

#include "FieldRecorder.h"

#include "../Utility/println.h"

#include <stdexcept>

// MARK:- P1D::FieldRecorder
//
std::string P1D::FieldRecorder::filepath(std::string const &wd, long const step_count) const
{
    constexpr char    prefix[] = "field";
    std::string const filename = std::string{prefix} + "-" + std::to_string(step_count) + ".csv";
    return is_master() ? wd + "/" + filename : null_dev;
}

P1D::FieldRecorder::FieldRecorder(parallel::mpi::Comm _comm)
: Recorder{Input::field_recording_frequency, std::move(_comm)}
{
    // configure output stream
    //
    os.setf(os.scientific);
    os.precision(15);
}

void P1D::FieldRecorder::record(const Domain &domain, const long step_count)
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
    print(os, "Nx = ", domain.params.Nx, '\n');
    //
    print(os, "dB1, dB2, dB3") << ", ";
    print(os, "dE1, dE2, dE3") << '\n';

    // contents
    //
    auto printer = [&os = this->os](Vector const &v) -> std::ostream & {
        return print(os, v.x, ", ", v.y, ", ", v.z);
    };

    if (is_master()) {
        auto const bfield = *comm.gather<1>({domain.bfield.begin(), domain.bfield.end()}, master);
        auto const efield = *comm.gather<1>({domain.efield.begin(), domain.efield.end()}, master);

        if (bfield.size() != domain.params.Nx)
            throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                     + " - incorrect gathered bfield size"};
        if (efield.size() != domain.params.Nx)
            throw std::runtime_error{std::string{__PRETTY_FUNCTION__}
                                     + " - incorrect gathered efield size"};

        for (unsigned long i = 0; i < domain.params.Nx; ++i) {
            printer(domain.geomtr.cart2fac(bfield[i] - domain.geomtr.B0)) << ", ";
            printer(domain.geomtr.cart2fac(efield[i])) << '\n';
        }
    } else {
        comm.gather<1>(domain.bfield.begin(), domain.bfield.end(), nullptr, master);
        comm.gather<1>(domain.efield.begin(), domain.efield.end(), nullptr, master);
    }

    os.close();
}
