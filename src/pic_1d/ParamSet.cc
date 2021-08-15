/*
 * Copyright (c) 2020-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "ParamSet.h"

#include <cmath>
#include <map>
#include <variant>

PIC1D_BEGIN_NAMESPACE
ParamSet::ParamSet(unsigned const rank, Options const &opts)
: geomtr{ Input::O0, Input::theta * M_PI / 180 }
{
    // domain extent
    //
    static_assert(Input::Nx % Input::number_of_subdomains == 0,
                  "Nx should be divisible by number_of_subdomains");
    Real const Mx = Input::Nx / Input::number_of_subdomains;
    domain_extent = { rank * Mx, Mx };

    // optional parameters
    //
    std::map<std::string_view, std::variant<long *, bool *, std::string *>> const map{
        { "wd", &working_directory },
        { "outer_Nt", &outer_Nt },
        { "save", &snapshot_save },
        { "load", &snapshot_load },
    };
    for (auto const &[key, val] : *opts) {
        std::visit(val, map.at(key));
    }
}

namespace {
template <class Object> decltype(auto) write_attr(Object &obj, ParamSet const &params)
{
    using hdf5::make_type;
    using hdf5::Space;
    { // global parameters
        obj.attribute("number_of_worker_threads", make_type(params.number_of_worker_threads),
                      Space::scalar())
            .write(params.number_of_worker_threads);
        obj.attribute("number_of_subdomains", make_type(params.number_of_subdomains),
                      Space::scalar())
            .write(params.number_of_subdomains);
        obj.attribute("number_of_particle_parallelism",
                      make_type(params.number_of_particle_parallelism), Space::scalar())
            .write(params.number_of_particle_parallelism);
        obj.attribute("is_electrostatic", make_type<int>(), Space::scalar())
            .template write<int>(params.is_electrostatic);
        obj.attribute("c", make_type(params.c), Space::scalar()).write(params.c);
        obj.attribute("O0", make_type(params.O0), Space::scalar()).write(params.O0);
        obj.attribute("theta", make_type(params.theta), Space::scalar()).write(params.theta);
        obj.attribute("Dx", make_type(params.Dx), Space::scalar()).write(params.Dx);
        obj.attribute("Nx", make_type(params.Nx), Space::scalar()).write(params.Nx);
        obj.attribute("dt", make_type(params.dt), Space::scalar()).write(params.dt);
        obj.attribute("innerNt", make_type(params.inner_Nt), Space::scalar())
            .write(params.inner_Nt);
        obj.attribute("partNs", make_type<long>(), Space::scalar())
            .template write<long>(std::tuple_size_v<decltype(params.part_descs)>);
        obj.attribute("coldNs", make_type<long>(), Space::scalar())
            .template write<long>(std::tuple_size_v<decltype(params.cold_descs)>);
    }
    return obj;
}
} // namespace
auto operator<<(hdf5::Group &obj, ParamSet const &params) -> decltype(obj)
{
    return write_attr(obj, params);
}
auto operator<<(hdf5::Dataset &obj, ParamSet const &params) -> decltype(obj)
{
    return write_attr(obj, params);
}
PIC1D_END_NAMESPACE
