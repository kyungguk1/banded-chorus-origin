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

#ifndef Driver_h
#define Driver_h

#include "./Boundary/MasterDelegate.h"
#include "./Boundary/SubdomainDelegate.h"
#include "./Core/Domain.h"
#include "./ParamSet.h"
#include "./Recorder/Recorder.h"

#include <array>
#include <future>
#include <map>
#include <memory>
#include <string>

HYBRID1D_BEGIN_NAMESPACE
class [[nodiscard]] Driver {
    long                                             iteration_count{};
    unsigned const                                   rank, size;
    ParamSet const                                   params;
    std::unique_ptr<Domain>                          domain;
    std::unique_ptr<MasterDelegate>                  master;
    std::unique_ptr<SubdomainDelegate>               delegate;
    std::map<std::string, std::unique_ptr<Recorder>> recorders;

    struct [[nodiscard]] Worker {
        long                    iteration_count;
        Driver const *          driver;
        WorkerDelegate *        delegate;
        std::future<void>       handle;
        std::unique_ptr<Domain> domain;
        //
        void operator()();
        Worker() noexcept          = default;
        Worker(Worker &&) noexcept = default;
    };
    std::array<Worker, ParamSet::number_of_particle_parallelism - 1> workers;

public:
    ~Driver();
    Driver(unsigned rank, unsigned size, ParamSet const &params);
    Driver(Driver &&) = default;

    void operator()();

private:
    void master_loop();

    [[nodiscard]] static std::unique_ptr<Domain> make_domain(ParamSet const &params,
                                                             Delegate *      delegate);
};
HYBRID1D_END_NAMESPACE

#endif /* Driver_h */
