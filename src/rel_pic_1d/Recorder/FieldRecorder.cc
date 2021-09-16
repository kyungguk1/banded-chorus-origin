/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "FieldRecorder.h"

#include <algorithm>
#include <stdexcept>

PIC1D_BEGIN_NAMESPACE
std::string FieldRecorder::filepath(std::string const &wd, long const step_count) const
{
    if (!is_master())
        throw std::domain_error{ __PRETTY_FUNCTION__ };

    constexpr char    prefix[] = "field";
    std::string const filename = std::string{ prefix } + "-" + std::to_string(step_count) + ".h5";
    return wd + "/" + filename;
}

FieldRecorder::FieldRecorder(parallel::mpi::Comm _comm)
: Recorder{ Input::field_recording_frequency, std::move(_comm) }
{
}

void FieldRecorder::record(const Domain &domain, const long step_count)
{
    if (step_count % recording_frequency)
        return;

    if (is_master())
        record_master(domain, step_count);
    else
        record_worker(domain, step_count);
}

template <class Object>
decltype(auto) FieldRecorder::write_attr(Object &&obj, Domain const &domain, long const step)
{
    obj << domain.params;
    obj.attribute("step", hdf5::make_type(step), hdf5::Space::scalar()).write(step);

    auto const time = step * domain.params.dt;
    obj.attribute("time", hdf5::make_type(time), hdf5::Space::scalar()).write(time);

    return std::forward<Object>(obj);
}
template <class T>
auto FieldRecorder::write_data(std::vector<T> payload, hdf5::Group &root, char const *name)
{
    auto const [mspace, fspace] = get_space(payload);
    auto const type             = hdf5::make_type<Real>();
    auto       dset             = root.dataset(name, type, fspace);
    dset.write(fspace, payload.data(), type, mspace);
    return dset;
}

void FieldRecorder::record_master(const Domain &domain, const long step_count)
{
    std::string const path = filepath(domain.params.working_directory, step_count);

    // create hdf file and root group
    auto root = hdf5::File(hdf5::File::trunc_tag{}, path.c_str())
                    .group("field", hdf5::PList::gapl(), hdf5::PList::gcpl());

    // attributes
    write_attr(root, domain, step_count);

    // datasets
    auto const writer = [](auto payload, auto &root, auto *name) {
        return write_data(std::move(payload), root, name);
    };

    if (auto obj = comm.gather<1>(cart2fac(domain.bfield, domain.params.geomtr), master)
                       .unpack(writer, root, "B"))
        write_attr(std::move(obj), domain, step_count) << domain.bfield;

    if (auto obj = comm.gather<1>(cart2fac(domain.efield, domain.params.geomtr), master)
                       .unpack(writer, root, "E"))
        write_attr(std::move(obj), domain, step_count) << domain.efield;

    root.flush();
}
void FieldRecorder::record_worker(const Domain &domain, const long)
{
    comm.gather<1>(cart2fac(domain.bfield, domain.params.geomtr), master).unpack([](auto) {});
    comm.gather<1>(cart2fac(domain.efield, domain.params.geomtr), master).unpack([](auto) {});
}

auto FieldRecorder::cart2fac(BField const &bfield, Geometry const &geomtr) -> std::vector<Vector>
{
    std::vector<Vector> dB(bfield.size());
    std::transform(bfield.begin(), bfield.end(), begin(dB), [&geomtr](auto const &v) {
        return geomtr.cart2fac(v - geomtr.B0);
    });

    return dB;
}
auto FieldRecorder::cart2fac(EField const &efield, Geometry const &geomtr) -> std::vector<Vector>
{
    std::vector<Vector> dE(efield.size());
    std::transform(efield.begin(), efield.end(), begin(dE), [&geomtr](auto const &v) {
        return geomtr.cart2fac(v);
    });

    return dE;
}
PIC1D_END_NAMESPACE