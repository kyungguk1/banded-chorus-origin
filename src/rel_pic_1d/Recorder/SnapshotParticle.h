/*
 * Copyright (c) 2020-2022, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include "../Core/Domain.h"
#include <PIC/TypeMaps.h>

#include <HDF5Kit/HDF5Kit.h>
#include <ParallelKit/ParallelKit.h>
#include <limits>
#include <string>
#include <string_view>

PIC1D_BEGIN_NAMESPACE
class SnapshotParticle {
    using Particle            = RelativisticParticle;
    using interprocess_comm_t = parallel::Communicator<Particle, long>;
    using rank_t              = parallel::mpi::Rank;

    interprocess_comm_t const comm;
    std::size_t const         signature;
    std::string_view const    wd; // working directory

    static constexpr unsigned chunk_size = std::numeric_limits<short>::max();
    static constexpr auto     tag        = parallel::mpi::Tag{ 5679 };
    static constexpr auto     tag2       = tag + 1;
    static constexpr rank_t   master{ 0 };

    [[nodiscard]] bool is_master() const { return master == comm->rank(); }
    [[nodiscard]] auto filepath() const;

public:
    SnapshotParticle(parallel::mpi::Comm distributed_comm, ParamSet const &);

    // load/save interface
    friend void save(SnapshotParticle &&snapshot, Domain const &domain, long step_count)
    {
        (snapshot.*snapshot.save)(domain, step_count);
        snapshot.comm.barrier();
    }
    [[nodiscard]] friend long load(SnapshotParticle &&snapshot, Domain &domain)
    {
        auto const result = (snapshot.*snapshot.load)(domain);
        return snapshot.comm.barrier(), result;
    }

private:
    void (SnapshotParticle::*save)(Domain const &, long) const &;
    long (SnapshotParticle::*load)(Domain &) const &;

    void save_helper(hdf5::Group &root, PartSpecies const &, std::string const &basename) const;
    void save_master(Domain const &, long step_count) const &;
    void save_worker(Domain const &, long step_count) const &;

    class LoaderPredicate;
    void               load_helper(hdf5::Group const &root, PartSpecies &, std::string const &basename) const;
    auto               distribute_particles(PartSpecies &, LoaderPredicate &) const -> unsigned long;
    [[nodiscard]] long load_master(Domain &) const &;
    [[nodiscard]] long load_worker(Domain &) const &;
};
PIC1D_END_NAMESPACE
