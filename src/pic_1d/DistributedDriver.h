/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include "Boundary/Delegate.h"
#include "Core/Domain.h"
#include "ParamSet.h"
#include "Recorder/Recorder.h"

#include <ParallelKit/ParallelKit.h>
#include <map>
#include <memory>
#include <string>

PIC1D_BEGIN_NAMESPACE
namespace Distributed {
class [[nodiscard]] Driver {
    long                      iteration_count{};
    ParamSet const            params;
    parallel::mpi::Comm const world;
    parallel::mpi::Comm       subdomain_comm;
    parallel::mpi::Comm       distributed_particle_comm;
    std::unique_ptr<Domain>   domain;
    std::unique_ptr<Delegate> subdomain_delegate;
    std::unique_ptr<Delegate> distributed_particle_delegate;
    // std::map<std::string, std::unique_ptr<Recorder>> recorders;

public:
    ~Driver();
    Driver(parallel::mpi::Comm comm, ParamSet const &params);
    Driver(Driver const &) = delete;

    void operator()();

private:
    void update();

    [[nodiscard]] static std::unique_ptr<Domain> make_domain(ParamSet const &params, Delegate *delegate);
};
} // namespace Distributed
PIC1D_END_NAMESPACE
