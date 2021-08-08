/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "FieldRecorder.h"

#include "../Utility/TypeMaps.h"

#include <algorithm>
#include <stdexcept>

// MARK:- P1D::FieldRecorder
//
std::string P1D::FieldRecorder::filepath(std::string const &wd, long const step_count) const
{
    if (!is_master())
        throw std::domain_error{ __PRETTY_FUNCTION__ };

    constexpr char    prefix[] = "field";
    std::string const filename = std::string{ prefix } + "-" + std::to_string(step_count) + ".h5";
    return wd + "/" + filename;
}

P1D::FieldRecorder::FieldRecorder(parallel::mpi::Comm _comm)
: Recorder{ Input::field_recording_frequency, std::move(_comm) }
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

template <class Object>
decltype(auto) P1D::FieldRecorder::write_attr(Object &&obj, Domain const &domain, long const step)
{
    obj << domain.params;
    obj.attribute("step", hdf5::make_type(step), hdf5::Space::scalar()).write(step);

    auto const time = step * domain.params.dt;
    obj.attribute("time", hdf5::make_type(time), hdf5::Space::scalar()).write(time);

    return std::forward<Object>(obj);
}
auto P1D::FieldRecorder::write_data(std::vector<Vector> payload, hdf5::Group &root,
                                    char const *name)
{
    auto space = hdf5::Space::simple(payload.size());
    auto dset  = root.dataset(name, hdf5::make_type<Vector>(), space);

    space.select_all();
    dset.write(space, payload.data(), space);

    return dset;
}

void P1D::FieldRecorder::record_master(const Domain &domain, const long step_count)
{
    std::string const path = filepath(domain.params.working_directory, step_count);

    // create hdf file and root group
    auto root = hdf5::File(hdf5::File::trunc_tag{}, path.c_str())
                    .group("field", hdf5::PList::gapl(), hdf5::PList::gcpl());

    // attributes
    write_attr(root, domain, step_count);

    // datasets
    if (auto obj = comm.gather<1>(cart2fac(domain.bfield, domain.geomtr), master)
                       .unpack(&write_data, root, "B"))
        write_attr(std::move(obj), domain, step_count) << domain.bfield;

    if (auto obj = comm.gather<1>(cart2fac(domain.efield, domain.geomtr), master)
                       .unpack(&write_data, root, "E"))
        write_attr(std::move(obj), domain, step_count) << domain.efield;

    root.flush();
}
void P1D::FieldRecorder::record_worker(const Domain &domain, const long)
{
    comm.gather<1>(cart2fac(domain.bfield, domain.geomtr), master).unpack([](auto) {});
    comm.gather<1>(cart2fac(domain.efield, domain.geomtr), master).unpack([](auto) {});
}

auto P1D::FieldRecorder::cart2fac(BField const &bfield, Geometry const &geomtr)
    -> std::vector<Vector>
{
    std::vector<Vector> dB(bfield.size());
    std::transform(bfield.begin(), bfield.end(), begin(dB), [&geomtr](auto const &v) {
        return geomtr.cart2fac(v - geomtr.B0);
    });

    return dB;
}
auto P1D::FieldRecorder::cart2fac(EField const &efield, Geometry const &geomtr)
    -> std::vector<Vector>
{
    std::vector<Vector> dE(efield.size());
    std::transform(efield.begin(), efield.end(), begin(dE), [&geomtr](auto const &v) {
        return geomtr.cart2fac(v);
    });

    return dE;
}
