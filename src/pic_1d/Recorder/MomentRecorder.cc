/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "MomentRecorder.h"

#include <algorithm>
#include <stdexcept>

PIC1D_BEGIN_NAMESPACE
std::string MomentRecorder::filepath(std::string const &wd, long const step_count) const
{
    if (!is_master())
        throw std::domain_error{ __PRETTY_FUNCTION__ };

    constexpr char    prefix[] = "moment";
    std::string const filename = std::string{ prefix } + "-" + std::to_string(step_count) + ".h5";
    return wd + "/" + filename;
}

MomentRecorder::MomentRecorder(parallel::mpi::Comm _comm)
: Recorder{ Input::moment_recording_frequency, std::move(_comm) }
{
}

void MomentRecorder::record(const Domain &domain, const long step_count)
{
    if (step_count % recording_frequency)
        return;

    if (is_master())
        record_master(domain, step_count);
    else
        record_worker(domain, step_count);
}

template <class Object>
decltype(auto) MomentRecorder::write_attr(Object &&obj, Domain const &domain, long const step)
{
    obj << domain.params;
    obj.attribute("step", hdf5::make_type(step), hdf5::Space::scalar()).write(step);

    auto const time = step * domain.params.dt;
    obj.attribute("time", hdf5::make_type(time), hdf5::Space::scalar()).write(time);

    return std::forward<Object>(obj);
}
auto MomentRecorder::get_space(std::vector<Scalar> const &payload)
{
    static_assert(sizeof(Scalar) % sizeof(Real) == 0);
    static_assert(sizeof(Scalar) / sizeof(Real) == 1);

    auto mspace = hdf5::Space::simple(payload.size());
    mspace.select_all();

    auto fspace = hdf5::Space::simple(payload.size());
    fspace.select_all();

    return std::make_pair(mspace, fspace);
}
auto MomentRecorder::get_space(std::vector<Vector> const &payload)
{
    constexpr auto size = 3U;
    static_assert(sizeof(Vector) % sizeof(Real) == 0);
    static_assert(sizeof(Vector) / sizeof(Real) >= size);

    auto mspace = hdf5::Space::simple({ payload.size(), sizeof(Vector) / sizeof(Real) });
    mspace.select(H5S_SELECT_SET, { 0U, 0U }, { payload.size(), size });

    auto fspace = hdf5::Space::simple({ payload.size(), size });
    fspace.select_all();

    return std::make_pair(mspace, fspace);
}
auto MomentRecorder::get_space(std::vector<Tensor> const &payload)
{
    static_assert(sizeof(Tensor) % sizeof(Real) == 0);
    static_assert(sizeof(Tensor) / sizeof(Real) == 8);

    auto mspace = hdf5::Space::simple({ payload.size(), sizeof(Tensor) / sizeof(Real) });
    // diagonal
    mspace.select(H5S_SELECT_SET, { 0U, 0U }, { payload.size(), 3U });
    // off-diag
    mspace.select(H5S_SELECT_OR, { 0U, 4U }, { payload.size(), 3U });

    auto fspace = hdf5::Space::simple({ payload.size(), 6U });
    fspace.select_all();

    return std::make_pair(mspace, fspace);
}
template <class T>
auto MomentRecorder::write_data(std::vector<T> payload, hdf5::Group &root, char const *name)
{
    auto const [mspace, fspace] = get_space(payload);
    auto const type             = hdf5::make_type<Real>();
    auto       dset             = root.dataset(name, type, fspace);
    dset.write(fspace, payload.data(), type, mspace);
    return dset;
}

void MomentRecorder::record_master(const Domain &domain, long const step_count)
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
    unsigned idx   = 0;
    auto     label = [&idx](std::string const &prefix) {
        return prefix + '_' + std::to_string(idx);
    };
    for (unsigned i = 0; i < part_Ns; ++i, ++idx) {
        PartSpecies const &sp = domain.part_species.at(i);

        if (auto obj = comm.gather<0>({ sp.moment<0>().begin(), sp.moment<0>().end() }, master)
                           .unpack(&write_data<Scalar>, root, label("n").c_str())) {
            write_attr(std::move(obj), domain, step_count) << sp;
        }
        if (auto obj = comm.gather<1>(cart2fac(sp.moment<1>(), domain.geomtr), master)
                           .unpack(&write_data<Vector>, root, label("nV").c_str())) {
            write_attr(std::move(obj), domain, step_count) << sp;
        }
        if (auto obj = comm.gather<1>(cart2fac(sp.moment<2>(), domain.geomtr), master)
                           .unpack(&write_data<Vector>, root, label("nvv").c_str())) {
            write_attr(std::move(obj), domain, step_count) << sp;
        }
    }

    for (unsigned i = 0; i < cold_Ns; ++i, ++idx) {
        ColdSpecies const &sp = domain.cold_species.at(i);

        if (auto obj = comm.gather<0>({ sp.moment<0>().begin(), sp.moment<0>().end() }, master)
                           .unpack(&write_data<Scalar>, root, label("n").c_str())) {
            write_attr(std::move(obj), domain, step_count) << sp;
        }
        if (auto obj = comm.gather<1>(cart2fac(sp.moment<1>(), domain.geomtr), master)
                           .unpack(&write_data<Vector>, root, label("nV").c_str())) {
            write_attr(std::move(obj), domain, step_count) << sp;
        }
        if (auto obj = comm.gather<1>(cart2fac(sp.moment<2>(), domain.geomtr), master)
                           .unpack(&write_data<Vector>, root, label("nvv").c_str())) {
            write_attr(std::move(obj), domain, step_count) << sp;
        }
    }

    root.flush();
}
void MomentRecorder::record_worker(const Domain &domain, long)
{
    for (PartSpecies const &sp : domain.part_species) {
        comm.gather<0>(sp.moment<0>().begin(), sp.moment<0>().end(), nullptr, master);
        comm.gather<1>(cart2fac(sp.moment<1>(), domain.geomtr), master).unpack([](auto) {});
        comm.gather<1>(cart2fac(sp.moment<2>(), domain.geomtr), master).unpack([](auto) {});
    }
    for (ColdSpecies const &sp : domain.cold_species) {
        comm.gather<0>(sp.moment<0>().begin(), sp.moment<0>().end(), nullptr, master);
        comm.gather<1>(cart2fac(sp.moment<1>(), domain.geomtr), master).unpack([](auto) {});
        comm.gather<1>(cart2fac(sp.moment<2>(), domain.geomtr), master).unpack([](auto) {});
    }
}

auto MomentRecorder::cart2fac(VectorGrid const &mom1, Geometry const &geomtr) -> std::vector<Vector>
{
    std::vector<Vector> nV(mom1.size());
    std::transform(mom1.begin(), mom1.end(), begin(nV), [&geomtr](auto const &v) {
        return geomtr.cart2fac(v);
    });

    return nV;
}
auto MomentRecorder::cart2fac(TensorGrid const &mom2, Geometry const &geomtr) -> std::vector<Vector>
{
    std::vector<Vector> nvv(mom2.size());
    std::transform(mom2.begin(), mom2.end(), begin(nvv), [&geomtr](auto const &v) {
        return geomtr.cart2fac(v).lo();
    });

    return nvv;
}
PIC1D_END_NAMESPACE
