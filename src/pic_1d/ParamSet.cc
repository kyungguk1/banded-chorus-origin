/*
 * Copyright (c) 2020-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "ParamSet.h"

#include <array>
#include <map>
#include <stdexcept>
#include <variant>

PIC1D_BEGIN_NAMESPACE
ParamSet::ParamSet(long const subdomain_rank, Options const &opts)
: geomtr{ Input::xi, Input::Dx, Input::O0 }
{
    if (subdomain_rank < 0 || subdomain_rank >= Input::number_of_subdomains)
        throw std::invalid_argument{ std::string{ __PRETTY_FUNCTION__ } + " - invalid rank" };

    constexpr auto full_to_half_grid_shift = 0.5;
    // whole domain extent
    full_grid_whole_domain_extent = { -0.5 * Input::Nx, Input::Nx };
    half_grid_whole_domain_extent = full_grid_whole_domain_extent + full_to_half_grid_shift;

    if (!geomtr.is_valid(CurviCoord{ full_grid_whole_domain_extent.min() })
        || !geomtr.is_valid(CurviCoord{ full_grid_whole_domain_extent.max() })
        || !geomtr.is_valid(CurviCoord{ half_grid_whole_domain_extent.min() })
        || !geomtr.is_valid(CurviCoord{ half_grid_whole_domain_extent.max() }))
        throw std::invalid_argument(std::string{ __PRETTY_FUNCTION__ } + " - invalid domain extent");

    // subdomain extent
    auto const Mx     = Input::Nx / Input::number_of_subdomains;
    auto const offset = subdomain_rank * Mx;

    full_grid_subdomain_extent = { full_grid_whole_domain_extent.min() + offset, Mx };
    half_grid_subdomain_extent = full_grid_subdomain_extent + full_to_half_grid_shift;

    // optional parameters
    //
    std::map<std::string_view, std::variant<long *, bool *, std::string *>> const map{
        { "wd", &working_directory },
        { "outer_Nt", &outer_Nt },
        { "save", &snapshot_save },
        { "load", &snapshot_load },
        { "record_particle_at_init", &record_particle_at_init },
    };
    for (auto const &[key, val] : *opts) {
        std::visit(val, map.at(key));
    }
}

namespace {
template <class Object>
decltype(auto) write_attr(Object &obj, ParamSet const &params)
{
    using hdf5::make_type;
    using hdf5::Space;
    { // global parameters
        std::array const extent{ params.full_grid_whole_domain_extent.min(), params.full_grid_whole_domain_extent.max() };
        auto const       type = make_type(extent.front());
        obj.attribute("full_grid_domain_extent", type, Space::simple(extent.size()))
            .write(extent.data(), type);
        obj.attribute("number_of_worker_threads", make_type(params.number_of_worker_threads), Space::scalar())
            .write(params.number_of_worker_threads);
        obj.attribute("number_of_subdomains", make_type(params.number_of_subdomains), Space::scalar())
            .write(params.number_of_subdomains);
        obj.attribute("number_of_distributed_particle_subdomain_clones", make_type(params.number_of_distributed_particle_subdomain_clones), Space::scalar())
            .write(params.number_of_distributed_particle_subdomain_clones);
        obj.attribute("number_of_mpi_processes", make_type(params.number_of_mpi_processes), Space::scalar())
            .write(params.number_of_mpi_processes);
        obj.attribute("number_of_particle_parallelism", make_type(params.number_of_particle_parallelism), Space::scalar())
            .write(params.number_of_particle_parallelism);
        obj.attribute("is_electrostatic", make_type<int>(), Space::scalar())
            .template write<int>(params.is_electrostatic);
        obj.attribute("c", make_type(params.c), Space::scalar()).write(params.c);
        obj.attribute("O0", make_type(params.O0), Space::scalar()).write(params.O0);
        obj.attribute("xi", make_type(params.xi), Space::scalar()).write(params.xi);
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
