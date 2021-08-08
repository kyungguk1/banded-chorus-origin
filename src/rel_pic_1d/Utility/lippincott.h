/*
 * Copyright (c) 2020-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#ifndef lippincott_h
#define lippincott_h

#include "../Macros.h"

#include <ParallelKit/ParallelKit.h>
#include <cstdlib>
#include <exception>
#include <string>

PIC1D_BEGIN_NAMESPACE
void print_backtrace();
PIC1D_END_NAMESPACE

// error handling routines
//
namespace {

[[noreturn, maybe_unused]] void fatal_error(char const *reason) noexcept
{
    std::puts(reason);
    P1D::print_backtrace();
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

#endif /* lippincott_h */
