/*
 * Copyright (c) 2020-2021, Kyungguk Min
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

#ifndef Snapshot_h
#define Snapshot_h

#include "../Core/Domain.h"
#include "../Utility/TypeMaps.h"

#include <HDF5Kit/HDF5Kit.h>
#include <ParallelKit/ParallelKit.h>
#include <memory>
#include <string>
#include <vector>

PIC1D_BEGIN_NAMESPACE
class Snapshot {
public:
    using interprocess_comm_t = parallel::Communicator<Scalar, Vector, Tensor, Particle, long>;
    using rank_t              = parallel::mpi::Rank;

    static constexpr parallel::mpi::Tag tag{599};

private:
    interprocess_comm_t const comm;
    std::size_t const         signature;
    std::string const         wd; // working directory

    static constexpr rank_t   master{0};
    [[nodiscard]] bool        is_master() const { return master == comm->rank(); }
    [[nodiscard]] std::string filepath() const;

public:
    explicit Snapshot(parallel::mpi::Comm comm, ParamSet const &params);

private: // load/save
    void (Snapshot::*save)(Domain const &domain, long step_count) const &;
    long (Snapshot::*load)(Domain &domain) const &;

    template <class T, long N>
    auto save_helper(hdf5::Group &root, GridQ<T, N> const &payload,
                     std::string const &basename) const -> hdf5::Dataset;
    void save_helper(hdf5::Group &root, PartSpecies const &payload) const;
    void save_master(Domain const &domain, long step_count) const &;
    void save_worker(Domain const &domain, long step_count) const &;

    template <class T, long N>
    auto               load_helper(hdf5::Group const &root, GridQ<T, N> &payload,
                                   std::string const &basename) const -> hdf5::Dataset;
    void               load_helper(hdf5::Group const &root, PartSpecies &payload) const;
    [[nodiscard]] long load_master(Domain &domain) const &;
    [[nodiscard]] long load_worker(Domain &domain) const &;

private: // load/save interface
    friend void save(Snapshot &&snapshot, Domain const &domain, long step_count)
    {
        return (snapshot.*snapshot.save)(domain, step_count);
    }
    [[nodiscard]] friend long load(Snapshot &&snapshot, Domain &domain)
    {
        return (snapshot.*snapshot.load)(domain);
    }
};
PIC1D_END_NAMESPACE

#endif /* Snapshot_h */
