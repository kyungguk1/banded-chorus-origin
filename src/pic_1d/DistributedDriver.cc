/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "DistributedDriver.h"
#include "Boundary/DistributedParticleDelegate.h"
#include "Boundary/SubdomainDelegate.h"
#include "Recorder/EnergyRecorder.h"
#include "Recorder/FieldRecorder.h"
#include "Recorder/MomentRecorder.h"
#include "Recorder/ParticleRecorder.h"
#include "Recorder/Snapshot.h"
#include "Recorder/VHistogramRecorder.h"
#include <PIC/lippincott.h>
#include <PIC/println.h>

#include <chrono>
#include <exception>
#include <iostream>

PIC1D_BEGIN_NAMESPACE
namespace Distributed {
namespace {
    template <class F, class... Args>
    [[nodiscard]] auto measure(F &&f, Args &&...args) -> std::chrono::duration<double>
    {
        static_assert(std::is_invocable_v<F &&, Args &&...>);
        auto const start = std::chrono::steady_clock::now();
        {
            std::invoke(std::forward<F>(f), std::forward<Args>(args)...);
        }
        auto const end = std::chrono::steady_clock::now();

        return end - start;
    }
} // namespace

Driver::~Driver()
{
}
Driver::Driver(parallel::mpi::Comm _comm, ParamSet const &params)
: params{ params }, world{ std::move(_comm) }
{
    try {
        if (!world)
            fatal_error(__PRETTY_FUNCTION__, " - invalid mpi::Comm object");

        if (auto const size = world.size(); size != params.number_of_mpi_processes)
            fatal_error(__PRETTY_FUNCTION__, " - the mpi world size (= ", std::to_string(size),
                        ") is not the same as params.number_of_mpi_processes (= ", std::to_string(params.number_of_mpi_processes), ')');

        auto const world_rank = world.rank();

        // group comm's
        // say, there are 3 subdomains
        // then the grouping of subdomain_comm's are {0, 1, 2}, {3, 4, 5}, ...
        // and the grouping of distributed_particle_comm's are {0, 3, 6, ...}, {1, 4, 7, ...}, and {2, 5, 8, ...}
        //
        subdomain_comm            = world.split(world_rank / long{ params.number_of_subdomains });
        distributed_particle_comm = world.split(world_rank % long{ params.number_of_subdomains });

        // TODO: Hopefully, this is a temporary restriction.
        if (1 != distributed_particle_comm.size())
            fatal_error(__PRETTY_FUNCTION__, " - distributed particle parallelization is currently disabled");

        // TODO: init recorders
        //
        // recorders["energy"]    = std::make_unique<EnergyRecorder>(comm.duplicated(), params);
        // recorders["fields"]    = std::make_unique<FieldRecorder>(comm.duplicated());
        // recorders["moment"]    = std::make_unique<MomentRecorder>(comm.duplicated());
        // recorders["vhists"]    = std::make_unique<VHistogramRecorder>(comm.duplicated());
        // recorders["particles"] = std::make_unique<ParticleRecorder>(comm.duplicated());

        // init delegates
        //
        subdomain_delegate            = std::make_unique<SubdomainDelegate>(subdomain_comm.duplicated());
        distributed_particle_delegate = std::make_unique<DistributedParticleDelegate>(distributed_particle_comm.duplicated(), subdomain_delegate.get());

        // init domain
        //
        if (0 == world_rank)
            println(std::cout, __FUNCTION__, "> initializing domain(s)");
        domain = make_domain(params, distributed_particle_delegate.get());

        // init particles or load snapshot
        //
        if (params.snapshot_load) {
            if (0 == world_rank)
                print(std::cout, "\tloading snapshots") << std::endl;
            // TODO: Construct Snapshot
            // iteration_count = load(Snapshot{ comm.duplicated(), params }, *domain);
        } else {
            if (0 == world_rank)
                print(std::cout, "\tinitializing particles") << std::endl;

            // TODO: Communicate to Species how to divide up
            for (PartSpecies &sp : domain->part_species) {
                sp.populate();
            }
            for (ColdSpecies &sp : domain->cold_species) {
                sp.populate();
            }

            /* TODO: Implement this
            if (params.record_particle_at_init) {
                // first, collect particle moments
                for (PartSpecies &sp : domain->part_species) {
                    sp.collect_all();
                    delegate->gather(*domain, sp);
                }

                // then, dump
                if (auto const &recorder = recorders.at("particles"))
                    recorder->record(*domain, iteration_count);
                if (auto const &recorder = recorders.at("vhists"))
                    recorder->record(*domain, iteration_count);
            }
             */
        }
    } catch (std::exception const &e) {
        fatal_error(__PRETTY_FUNCTION__, " :: ", e.what());
    }
}
auto Driver::make_domain(ParamSet const &params, Delegate *delegate) -> std::unique_ptr<Domain>
{
    return std::make_unique<Domain>(params, delegate);
}

void Driver::operator()()
try {
    auto const elapsed = measure(&Driver::update, this);
    if (0 == world.rank())
        println(std::cout, "%% time elapsed: ", elapsed.count(), 's');

    // TODO: take snapshot
    /*
    if (params.snapshot_save) {
        if (0 == world.rank())
            print(std::cout, "\tsaving snapshots") << std::endl;
        save(Snapshot{ world.duplicated(), params }, *domain, iteration_count);
    }
     */
} catch (std::exception const &e) {
    fatal_error(__PRETTY_FUNCTION__, " :: ", e.what());
}
void Driver::update()
try {
    for (long outer_step = 1; outer_step <= params.outer_Nt; ++outer_step) {
        if (0 == world.rank())
            println(std::cout, __FUNCTION__, "> ",
                    "steps(x", params.inner_Nt, ") = ", outer_step, "/", params.outer_Nt,
                    "; time = ", iteration_count * params.dt);

        // inner loop
        //
        domain->advance_by(params.inner_Nt);

        // increment step count
        //
        iteration_count += params.inner_Nt;

        // TODO: record data
        /*
        for (auto &pair : recorders) {
            if (pair.second)
                pair.second->record(*domain, iteration_count);
        }
         */
    }
} catch (std::exception const &e) {
    fatal_error(__PRETTY_FUNCTION__, " :: ", e.what());
}
} // namespace Distributed
PIC1D_END_NAMESPACE
