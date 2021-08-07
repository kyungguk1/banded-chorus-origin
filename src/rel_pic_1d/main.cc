/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "./Driver.h"
#include "./Utility/lippincott.h"
#include "./Utility/println.h"

#include <future>
#include <iostream>
#include <stdexcept>

int main(int argc, char *argv[])
try {
    using parallel::mpi::Comm;
    {
        constexpr bool enable_mpi_thread = false;
        int const      required = enable_mpi_thread ? MPI_THREAD_MULTIPLE : MPI_THREAD_SINGLE;
        int            provided;
        if (Comm::init(&argc, &argv, required, provided) != MPI_SUCCESS) {
            println(std::cout, "%% ", __PRETTY_FUNCTION__,
                    " - mpi::Comm::init(...) returned error");
            return 1;
        }
        if (provided < required) {
            println(std::cout, "%% ", __PRETTY_FUNCTION__, " - provided < required");
            return 1;
        }
    }

    if (auto world = Comm::world().duplicated()) {
        auto const rank = world.rank();
        auto const opts = P1D::Options{ { argv, argv + argc } };
        P1D::Driver{ std::move(world), { static_cast<unsigned>(rank), opts } }();
    } else {
        throw std::runtime_error{ std::string{ __PRETTY_FUNCTION__ } + " - invalid mpi::Comm" };
    }

    Comm::deinit();
    return 0;
} catch (...) {
    lippincott();
}
