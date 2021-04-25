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

#include "../Utility/TypeMaps.h"

#include <HDF5Kit/HDF5Kit.h>
#include <algorithm>
#include <iterator>
#include <stdexcept>

// MARK:- P1D::FieldRecorder
//
std::string P1D::FieldRecorder::filepath(std::string const &wd, long const step_count) const
{
    if (!is_master())
        throw std::domain_error{__PRETTY_FUNCTION__};

    constexpr char    prefix[] = "field";
    std::string const filename = std::string{prefix} + "-" + std::to_string(step_count) + ".h5";
    return wd + "/" + filename;
}

P1D::FieldRecorder::FieldRecorder(parallel::mpi::Comm _comm)
: Recorder{Input::field_recording_frequency, std::move(_comm)}
{
}

void P1D::FieldRecorder::record(const Domain &domain, const long step_count)
{
    if (step_count % recording_frequency)
        return;

    if (is_master())
        record_master(domain, step_count);
    else
        record_worker(domain, step_count);
}

namespace {
template <class Object>
decltype(auto) write_attr(Object &&obj, P1D::Domain const &domain, long const step)
{
    obj << domain.params;
    obj.attribute("step", hdf5::make_type(step), hdf5::Space::scalar()).write(step);

    auto const time = step * domain.params.dt;
    obj.attribute("time", hdf5::make_type(time), hdf5::Space::scalar()).write(time);

    return std::forward<Object>(obj);
}
template <class T> auto write_data(std::vector<T> payload, hdf5::Group &root, char const *name)
{
    auto space = hdf5::Space::simple(payload.size());
    auto dset  = root.dataset(name, hdf5::make_type<T>(), space);

    space.select_all();
    dset.write(space, payload.data(), space);

    return dset;
}
auto write_vector(std::vector<P1D::Vector> payload, P1D::Geometry const &geomtr, hdf5::Group &root,
                  char const *name)
{
    std::transform(begin(payload), end(payload), begin(payload), [&geomtr](auto const &v) {
        return geomtr.cart2fac(v);
    });
    return write_data(std::move(payload), root, name);
}
} // namespace
void P1D::FieldRecorder::record_master(const Domain &domain, const long step_count)
{
    std::string const path = filepath(domain.params.working_directory, step_count);

    // create hdf file and root group
    auto root = hdf5::File(hdf5::File::trunc_tag{}, path.c_str())
                    .group("field", hdf5::PList::gapl(), hdf5::PList::gcpl());

    // attributes
    write_attr(root, domain, step_count);

    // datasets
    std::vector<Vector> bfield{domain.bfield.begin(), domain.bfield.end()};
    std::transform(begin(bfield), end(bfield), begin(bfield),
                   [B0 = domain.geomtr.B0](auto const &v) {
                       return v - B0;
                   });
    if (auto obj = comm.gather<1>(bfield, master).unpack(&write_vector, domain.geomtr, root, "B"))
        write_attr(std::move(obj), domain, step_count) << domain.bfield;

    if (auto obj = comm.gather<1>({domain.efield.begin(), domain.efield.end()}, master)
                       .unpack(&write_vector, domain.geomtr, root, "E"))
        write_attr(std::move(obj), domain, step_count) << domain.efield;

    root.flush();
}
void P1D::FieldRecorder::record_worker(const Domain &domain, const long)
{
    std::vector<Vector> bfield{domain.bfield.begin(), domain.bfield.end()};
    std::transform(begin(bfield), end(bfield), begin(bfield),
                   [B0 = domain.geomtr.B0](auto const &v) {
                       return v - B0;
                   });
    comm.gather<1>(bfield.data(), std::next(bfield.data(), domain.bfield.size()), nullptr, master);
    comm.gather<1>(domain.efield.begin(), domain.efield.end(), nullptr, master);
}
