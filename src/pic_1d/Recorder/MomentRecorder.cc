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

#include "../Utility/TypeMaps.h"

#include <HDF5Kit/HDF5Kit.h>
#include <algorithm>
#include <stdexcept>

// MARK:- P1D::MomentRecorder
//
std::string P1D::MomentRecorder::filepath(std::string const &wd, long const step_count) const
{
    if (!is_master())
        throw std::domain_error{__PRETTY_FUNCTION__};

    constexpr char    prefix[] = "moment";
    std::string const filename = std::string{prefix} + "-" + std::to_string(step_count) + ".h5";
    return wd + "/" + filename;
}

P1D::MomentRecorder::MomentRecorder(parallel::mpi::Comm _comm)
: Recorder{Input::moment_recording_frequency, std::move(_comm)}
{
}

void P1D::MomentRecorder::record(const Domain &domain, const long step_count)
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
auto write_scalar(std::vector<P1D::Scalar> payload, hdf5::Group &root, char const *name)
{
    return write_data(std::move(payload), root, name);
}
auto write_vector(std::vector<P1D::Vector> payload, P1D::Geometry const &geomtr, hdf5::Group &root,
                  char const *name)
{
    std::transform(begin(payload), end(payload), begin(payload), [&geomtr](auto const &v) {
        return geomtr.cart2fac(v);
    });
    return write_data(std::move(payload), root, name);
}
auto write_tensor(std::vector<P1D::Tensor> payload, P1D::Geometry const &geomtr, hdf5::Group &root,
                  char const *name)
{
    std::vector<P1D::Vector> transformed(payload.size());
    std::transform(begin(payload), end(payload), begin(transformed), [&geomtr](auto const &v) {
        return geomtr.cart2fac(v);
    });
    return write_data(std::move(transformed), root, name);
}
} // namespace
void P1D::MomentRecorder::record_master(const Domain &domain, long const step_count)
{
    std::string const path = filepath(domain.params.working_directory, step_count);

    // create hdf file and root group
    auto root = hdf5::File(hdf5::File::trunc_tag{}, path.c_str())
                    .group("moment", hdf5::PList::gapl(), hdf5::PList::gcpl());

    // attributes
    auto const part_Ns = domain.part_species.size();
    auto const cold_Ns = domain.cold_species.size();
    auto const Ns      = part_Ns + cold_Ns;
    write_attr(root, domain, step_count)
        .attribute("Ns", hdf5::make_type(Ns), hdf5::Space::scalar())
        .write(Ns);

    // datasets
    unsigned idx = 0;
    auto label = [&idx](std::string const &prefix) {
        return prefix + "_" + std::to_string(idx);
    };
    for (unsigned i = 0; i < part_Ns; ++i, ++idx) {
        PartSpecies const &sp = domain.part_species.at(i);

        if (auto obj = comm.gather<0>({sp.moment<0>().begin(), sp.moment<0>().end()}, master)
                           .unpack(&write_scalar, root, label("n").c_str())) {
            write_attr(std::move(obj), domain, step_count) << sp;
        }
        if (auto obj = comm.gather<1>({sp.moment<1>().begin(), sp.moment<1>().end()}, master)
                           .unpack(&write_vector, domain.geomtr, root, label("nV").c_str())) {
            write_attr(std::move(obj), domain, step_count) << sp;
        }
        if (auto obj = comm.gather<2>({sp.moment<2>().begin(), sp.moment<2>().end()}, master)
                           .unpack(&write_tensor, domain.geomtr, root, label("nvv").c_str())) {
            write_attr(std::move(obj), domain, step_count) << sp;
        }
    }

    for (unsigned i = 0; i < cold_Ns; ++i, ++idx) {
        ColdSpecies const &sp = domain.cold_species.at(i);

        if (auto obj = comm.gather<0>({sp.moment<0>().begin(), sp.moment<0>().end()}, master)
                           .unpack(&write_scalar, root, label("n").c_str())) {
            write_attr(std::move(obj), domain, step_count) << sp;
        }
        if (auto obj = comm.gather<1>({sp.moment<1>().begin(), sp.moment<1>().end()}, master)
                           .unpack(&write_vector, domain.geomtr, root, label("nV").c_str())) {
            write_attr(std::move(obj), domain, step_count) << sp;
        }
        if (auto obj = comm.gather<2>({sp.moment<2>().begin(), sp.moment<2>().end()}, master)
                           .unpack(&write_tensor, domain.geomtr, root, label("nvv").c_str())) {
            write_attr(std::move(obj), domain, step_count) << sp;
        }
    }

    root.flush();
}
void P1D::MomentRecorder::record_worker(const Domain &domain, long const)
{
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
