/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "Recorder.h"

#include <limits>
#include <stdexcept>

namespace {
constexpr long large_int = std::numeric_limits<unsigned>::max();
}

// MARK:- P1D::Recorder
//
P1D::Recorder::Recorder(unsigned const recording_frequency, parallel::mpi::Comm _comm)
: recording_frequency{ recording_frequency ? recording_frequency * Input::inner_Nt : large_int }
, comm{ std::move(_comm) }
{
    if (!comm->operator bool())
        throw std::invalid_argument{ std::string{ __PRETTY_FUNCTION__ } + " - invalid mpi::Comm" };
}
