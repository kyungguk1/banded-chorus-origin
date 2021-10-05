/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "Driver.h"
#include "Recorder/EnergyRecorder.h"
#include "Recorder/FieldRecorder.h"
#include "Recorder/MomentRecorder.h"
#include "Recorder/ParticleRecorder.h"
#include "Recorder/Snapshot.h"
#include "Recorder/VHistogramRecorder.h"
#include <PIC/lippincott.h>
#include <PIC/println.h>

#include <chrono>
#include <iostream>

PIC1D_BEGIN_NAMESPACE
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
: comm{ std::move(_comm) }, params{ params }
{
    try {
        auto const &comm = this->comm;
        if (!comm)
            throw std::invalid_argument{ std::string{ __PRETTY_FUNCTION__ } + " - invalid mpi::Comm object" };

        if (auto const size = comm.size(); size != params.number_of_subdomains)
            throw std::runtime_error{ std::string{ __PRETTY_FUNCTION__ } + " - the mpi comm size (= " + std::to_string(size) + ") is not the same as number_of_subdomains (= " + std::to_string(params.number_of_subdomains) + ')' };

        auto const rank = comm.rank();

        // init recorders
        //
        recorders["energy"]    = std::make_unique<EnergyRecorder>(comm.duplicated(), params);
        recorders["fields"]    = std::make_unique<FieldRecorder>(comm.duplicated());
        recorders["moment"]    = std::make_unique<MomentRecorder>(comm.duplicated());
        recorders["vhists"]    = std::make_unique<VHistogramRecorder>(comm.duplicated());
        recorders["particles"] = std::make_unique<ParticleRecorder>(comm.duplicated());

        // init master delegate
        //
        delegate = std::make_unique<SubdomainDelegate>(comm.duplicated());
        master   = std::make_unique<MasterDelegate>(delegate.get());

        // init master domain
        //
        if (0 == rank)
            println(std::cout, __FUNCTION__, "> initializing domain(s)");
        domain = make_domain(params, master.get());

        // init particles or load snapshot
        //
        if (params.snapshot_load) {
            if (0 == rank)
                print(std::cout, "\tloading snapshots") << std::endl;
            iteration_count = load(Snapshot{ comm.duplicated(), params }, *domain);
        } else {
            if (0 == rank)
                print(std::cout, "\tinitializing particles") << std::endl;

            for (PartSpecies &sp : domain->part_species) {
                sp.populate();
            }
            for (ColdSpecies &sp : domain->cold_species) {
                sp.populate();
            }

            if (params.record_particle_at_init) {
                if (auto const &recorder = recorders.at("particles"))
                    recorder->record(*domain, iteration_count);
                if (auto const &recorder = recorders.at("vhists"))
                    recorder->record(*domain, iteration_count);
            }
        }
    } catch (...) {
        lippincott();
    }
}
auto Driver::make_domain(ParamSet const &params, Delegate *delegate) -> std::unique_ptr<Domain>
{
    return std::make_unique<Domain>(params, delegate);
}

void Driver::operator()()
try {
    // worker setup
    //
    for (unsigned i = 0; i < workers.size(); ++i) {
        Worker &worker         = workers[i];
        worker.driver          = this;
        worker.iteration_count = iteration_count;
        worker.delegate        = &master->workers.at(i);
        worker.domain          = make_domain(params, worker.delegate);
        worker.handle          = std::async(std::launch::async, worker.delegate->wrap_loop(std::ref(worker)), worker.domain.get());
    }

    // master loop
    //
    auto const elapsed = measure(master->wrap_loop(&Driver::master_loop, this), this->domain.get());
    if (0 == comm.rank())
        println(std::cout, "%% time elapsed: ", elapsed.count(), 's');

    // worker teardown
    //
    for (Worker &worker : workers) {
        worker.handle.get();
        worker.domain.reset();
    }

    // take snapshot
    //
    if (params.snapshot_save) {
        if (0 == comm.rank())
            print(std::cout, "\tsaving snapshots") << std::endl;
        save(Snapshot{ comm.duplicated(), params }, *domain, iteration_count);
    }
} catch (...) {
    lippincott();
}
void Driver::master_loop()
try {
    for (long outer_step = 1; outer_step <= domain->params.outer_Nt; ++outer_step) {
        if (0 == comm.rank())
            println(std::cout, __FUNCTION__, "> ",
                    "steps(x", domain->params.inner_Nt, ") = ", outer_step, "/", domain->params.outer_Nt,
                    "; time = ", iteration_count * domain->params.dt);

        // inner loop
        //
        domain->advance_by(domain->params.inner_Nt);

        // increment step count
        //
        iteration_count += domain->params.inner_Nt;

        // record data
        //
        if (iteration_count % this->recorders.at("vhists")->recording_frequency
            && iteration_count % this->recorders.at("particles")->recording_frequency) {
            // no particle collection needed
            //
            for (auto &pair : recorders) {
                if (pair.second)
                    pair.second->record(*domain, iteration_count);
            }
        } else {
            // collect particles before recording
            //
            auto const *delegate = master.get();
            delegate->teardown(*domain);

            // record data
            //
            for (auto &pair : recorders) {
                if (pair.second)
                    pair.second->record(*domain, iteration_count);
            }

            // re-distribute particles
            //
            delegate->setup(*domain);
        }
    }
} catch (...) {
    lippincott();
}
void Driver::Worker::operator()()
try {
    for (long outer_step = 1; outer_step <= domain->params.outer_Nt; ++outer_step) {
        // inner loop
        //
        domain->advance_by(domain->params.inner_Nt);

        // increment step count
        //
        iteration_count += domain->params.inner_Nt;

        // record data
        //
        if (iteration_count % driver->recorders.at("vhists")->recording_frequency
            && iteration_count % driver->recorders.at("particles")->recording_frequency) {
            // no particle collection needed
            //
        } else {
            // collect particles before recording
            //
            delegate->teardown(*domain);
            delegate->setup(*domain);
        }
    }
} catch (...) {
    lippincott();
}
PIC1D_END_NAMESPACE
