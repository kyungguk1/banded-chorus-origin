/*
 * Copyright (c) 2020-2021, Kyungguk Min
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

#include "ParamSet.h"

#include <limits>
#include <map>
#include <string_view>
#include <variant>

namespace {
constexpr auto quiet_nan = std::numeric_limits<double>::quiet_NaN();
}
P1D::ParamSet::ParamSet() noexcept : domain_extent{ quiet_nan, quiet_nan }
{
}
P1D::ParamSet::ParamSet(unsigned const rank, Options const &opts) : ParamSet{}
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
template <class Object> decltype(auto) write_attr(Object &obj, P1D::ParamSet const &params)
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
auto P1D::operator<<(hdf5::Group &obj, P1D::ParamSet const &params) -> decltype(obj)
{
    return write_attr(obj, params);
}
auto P1D::operator<<(hdf5::Dataset &obj, P1D::ParamSet const &params) -> decltype(obj)
{
    return write_attr(obj, params);
}
