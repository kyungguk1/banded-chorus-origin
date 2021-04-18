//
// Copyright (c) 2020, Kyungguk Min
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// - Redistributions of source code must retain the above copyright notice, this
//  list of conditions and the following disclaimer.
//
// - Redistributions in binary form must reproduce the above copyright notice,
//  this list of conditions and the following disclaimer in the documentation
//  and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include "ParamSet.h"

#include <limits>
#include <map>
#include <string_view>
#include <variant>

namespace {
constexpr auto quiet_nan = std::numeric_limits<double>::quiet_NaN();
}
P1D::ParamSet::ParamSet() noexcept : domain_extent{quiet_nan, quiet_nan}
{
}
P1D::ParamSet::ParamSet(unsigned const rank, Options const &opts) : ParamSet{}
{
    // domain extent
    //
    static_assert(Input::Nx % Input::number_of_subdomains == 0,
                  "Nx should be divisible by number_of_subdomains");
    Real const Mx = Input::Nx / Input::number_of_subdomains;
    domain_extent = {rank * Mx, Mx};

    // optional parameters
    //
    std::map<std::string_view, std::variant<long *, bool *, std::string *>> const map{
        {"wd", &working_directory},
        {"outer_Nt", &outer_Nt},
        {"save", &snapshot_save},
        {"load", &snapshot_load}};
    for (auto const &[key, val] : *opts) {
        std::visit(val, map.at(key));
    }
}
