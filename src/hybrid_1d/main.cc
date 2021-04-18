/*
 * Copyright (c) 2019, Kyungguk Min
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

#include <array>
#include <functional>
#include <future>

#if defined(DEBUG)
#include "./Utility/Options.h"
#include "./VDF/BitReversedPattern.h"
#include "./VDF/LossconeVDF.h"
#endif

int main([[maybe_unused]] int argc, [[maybe_unused]] const char *argv[])
try {
    using namespace H1D;
    {
        constexpr unsigned size = Input::number_of_subdomains;
        auto               task = [opts = Options{{argv, argv + argc}}](unsigned const rank) {
            // construction of Driver should be done on their own thread
            return Driver{rank, size, {rank, opts}}();
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
    //
#if defined(DEBUG)
//    test_BitReversedPattern();
//    test_option_parser();
//    test_LossconeVDF();
#endif
    //
    return 0;
} catch (...) {
    lippincott();
}
