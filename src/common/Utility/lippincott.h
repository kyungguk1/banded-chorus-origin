/*
 * Copyright (c) 2020-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#ifndef COMMON_LIPPINCOTT_h
#define COMMON_LIPPINCOTT_h

#include <common-config.h>

#include <ParallelKit/ParallelKit.h>
#include <cstdlib>
#include <exception>
#include <string>

// error handling routines
//
COMMON_BEGIN_NAMESPACE
void print_backtrace();
COMMON_END_NAMESPACE

namespace {

[[noreturn, maybe_unused]] void fatal_error(char const *reason) noexcept
{
    std::puts(reason);
    COMMON_NAMESPACE::print_backtrace();
    if (parallel::mpi::Comm::is_initialized())
        MPI_Abort(MPI_COMM_WORLD, 1);
    std::abort();
}
[[noreturn, maybe_unused]] void fatal_error(std::string const &reason) noexcept
{
    fatal_error(reason.c_str());
}
[[noreturn, maybe_unused]] void lippincott() noexcept
try {
    throw;
} catch (std::exception const &e) {
    fatal_error(e.what());
} catch (...) {
    fatal_error("Unknown exception");
}

} // namespace

#endif /* COMMON_LIPPINCOTT_h */
