/*
 * Copyright (c) 2020, Kyungguk Min
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

#ifndef ParamSet_h
#define ParamSet_h

#include "./InputWrapper.h"
#include "./Utility/Options.h"

HYBRID1D_BEGIN_NAMESPACE
struct [[nodiscard]] ParamSet : public Input {

    /// number of threads for particle async update
    ///
    static constexpr unsigned number_of_particle_parallelism
        = (number_of_worker_threads + 1) / number_of_subdomains;

    /// index sequence of kinetic plasma descriptors
    ///
    using part_indices = std::make_index_sequence<std::tuple_size_v<decltype(part_descs)>>;

    /// index sequence of cold plasma descriptors
    ///
    using cold_indices = std::make_index_sequence<std::tuple_size_v<decltype(cold_descs)>>;

public:
    Range       domain_extent;
    long        outer_Nt{Input::outer_Nt};
    std::string working_directory{Input::working_directory};
    bool        snapshot_save{false};
    bool        snapshot_load{false};
    //
    ParamSet() noexcept;
    ParamSet(unsigned rank, Options const &opts);

private:
    template <class... Ts, class Int, Int... Is>
    [[nodiscard]] static constexpr auto _serialize(std::tuple<Ts...> const &t,
                                                   std::integer_sequence<Int, Is...>) noexcept
    {
        return std::tuple_cat(serialize(std::get<Is>(t))...);
    }
    [[nodiscard]] friend constexpr auto serialize(ParamSet const &params) noexcept
    {
        auto const global = std::make_tuple(params.algorithm, params.c, params.O0, params.theta,
                                            params.Dx, params.Nx, params.dt, params.inner_Nt);
        auto const efluid = serialize(params.efluid_desc);
        auto const parts  = _serialize(params.part_descs, part_indices{});
        auto const colds  = _serialize(params.cold_descs, cold_indices{});
        return std::tuple_cat(global, efluid, parts, colds);
    }
};
HYBRID1D_END_NAMESPACE

#endif /* ParamSet_h */
