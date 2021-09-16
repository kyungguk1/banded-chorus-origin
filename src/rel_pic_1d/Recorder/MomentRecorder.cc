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
    auto const writer = [](auto payload, auto &root, auto *name) {
        return write_data(std::move(payload), root, name);
    };

    unsigned idx = 0;
    for (unsigned i = 0; i < part_Ns; ++i, ++idx) {
        PartSpecies const &sp = domain.part_species.at(i);

        auto parent = [&root, name = std::to_string(idx)] {
            return root.group(name.c_str(), hdf5::PList::gapl(), hdf5::PList::gcpl());
        }();
        write_attr(parent, domain, step_count) << sp;

        comm.gather<0>({ sp.moment<0>().begin(), sp.moment<0>().end() }, master)
            .unpack(writer, parent, "n");
        comm.gather<1>(cart2fac(sp.moment<1>(), domain.params.geomtr), master)
            .unpack(writer, parent, "nV");
        comm.gather<2>(cart2fac(sp.moment<2>(), domain.params.geomtr), master)
            .unpack(writer, parent, "Mij");
    }
    for (unsigned i = 0; i < cold_Ns; ++i, ++idx) {
        ColdSpecies const &sp = domain.cold_species.at(i);

        auto parent = [&root, name = std::to_string(idx)] {
            return root.group(name.c_str(), hdf5::PList::gapl(), hdf5::PList::gcpl());
        }();
        write_attr(parent, domain, step_count) << sp;

        comm.gather<0>({ sp.moment<0>().begin(), sp.moment<0>().end() }, master)
            .unpack(writer, parent, "n");
        comm.gather<1>(cart2fac(sp.moment<1>(), domain.params.geomtr), master)
            .unpack(writer, parent, "nV");
        comm.gather<2>(cart2fac(sp.moment<2>(), domain.params.geomtr), master)
            .unpack(writer, parent, "Mij");
    }

    root.flush();
}
void MomentRecorder::record_worker(const Domain &domain, long)
{
    for (PartSpecies const &sp : domain.part_species) {
        comm.gather<0>(sp.moment<0>().begin(), sp.moment<0>().end(), nullptr, master);
        comm.gather<1>(cart2fac(sp.moment<1>(), domain.params.geomtr), master)
            .unpack([](auto) {});
        comm.gather<2>(cart2fac(sp.moment<2>(), domain.params.geomtr), master)
            .unpack([](auto) {});
    }
    for (ColdSpecies const &sp : domain.cold_species) {
        comm.gather<0>(sp.moment<0>().begin(), sp.moment<0>().end(), nullptr, master);
        comm.gather<1>(cart2fac(sp.moment<1>(), domain.params.geomtr), master)
            .unpack([](auto) {});
        comm.gather<2>(cart2fac(sp.moment<2>(), domain.params.geomtr), master)
            .unpack([](auto) {});
    }
}

auto MomentRecorder::cart2fac(VectorGrid const &mom1, Geometry const &geomtr) -> std::vector<Vector>
{
    std::vector<Vector> nV(mom1.size());
    std::transform(mom1.begin(), mom1.end(), begin(nV), [&geomtr](auto const &nV) {
        return geomtr.cart2fac(nV);
    });

    return nV;
}
auto MomentRecorder::cart2fac(FourTensorGrid const &mom2, Geometry const &geomtr) -> std::vector<FourTensor>
{
    std::vector<FourTensor> Mij(mom2.size());
    std::transform(mom2.begin(), mom2.end(), begin(Mij), [&geomtr](auto const &Mij) {
        return geomtr.cart2fac(Mij);
    });

    return Mij;
}
PIC1D_END_NAMESPACE