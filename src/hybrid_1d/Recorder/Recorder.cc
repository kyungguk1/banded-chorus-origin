//
// Copyright (c) 2019, Kyungguk Min
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

#include "Recorder.h"

#include <limits>
#include <stdexcept>

namespace {
constexpr long large_int = std::numeric_limits<unsigned>::max();
}

H1D::Recorder::message_dispatch_t H1D::Recorder::dispatch{H1D::ParamSet::number_of_subdomains};
H1D::Recorder::Recorder(unsigned const recording_frequency, unsigned const rank,
                        unsigned const size)
: recording_frequency{recording_frequency ? recording_frequency * Input::inner_Nt : large_int}
, comm{dispatch.comm(rank)}
, size{size}
, all_ranks{}
, all_but_master{}
{
    if (size > ParamSet::number_of_subdomains)
        throw std::invalid_argument{__PRETTY_FUNCTION__};

    for (unsigned rank = 0; is_master() && rank < size; ++rank) {
        all_ranks.emplace_back(rank);
        if (master != rank)
            all_but_master.emplace_back(rank);
    }
}
