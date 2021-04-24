/*
 * Copyright (c) 2019-2021, Kyungguk Min
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

#include "./Driver.h"
#include "./Utility/lippincott.h"
#include "./Utility/println.h"

#include <array>
#include <functional>
#include <future>
#include <iostream>

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
        } else if (provided < required) {
            println(std::cout, "%% ", __PRETTY_FUNCTION__, " - provided < required");
            return 1;
        }
    }

    using namespace P1D;
    if (Comm::world().size() > 1) { // mpi version for domain decomposition

        auto const opts  = Options{{argv, argv + argc}};
        auto       world = Comm::world().duplicated();
        auto const rank  = world.rank();
        mpi::Driver{std::move(world), {static_cast<unsigned>(rank), opts}}();

    } else { // multi-threaded version for domain decomposition
        constexpr unsigned size = Input::number_of_subdomains;
        auto               task = [opts = Options{{argv, argv + argc}}](unsigned const rank) {
            // construction of Driver should be done on their own thread
            return thread::Driver{rank, size, {rank, opts}}();
        };
        //
        std::array<std::future<void>, size> workers;
        std::packaged_task<void(unsigned)>  main_task{task};
        workers.at(0) = main_task.get_future();
        for (unsigned rank = 1; rank < size; ++rank) {
            workers.at(rank) = std::async(std::launch::async, task, rank);
        }
        std::invoke(main_task, 0);
        //
        for (auto &f : workers) {
            f.get();
        }
    }

    Comm::deinit();
    return 0;
} catch (...) {
    lippincott();
}
