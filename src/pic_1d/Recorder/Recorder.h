/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include "../Core/Domain.h"
#include <PIC/TypeMaps.h>

#include <HDF5Kit/HDF5Kit.h>
#include <ParallelKit/ParallelKit.h>
#include <type_traits>
#include <utility>
#include <vector>

PIC1D_BEGIN_NAMESPACE
class Recorder {
public:
    long const recording_frequency;

    virtual ~Recorder() = default;

    virtual void record(Domain const &domain, long step_count) = 0;

protected:
    Recorder(unsigned recording_frequency, parallel::mpi::Comm subdomain_comm, parallel::mpi::Comm distributed_particle_comm);

    using rank_t = parallel::mpi::Rank;

    parallel::Communicator<Scalar, Vector, Tensor> const subdomain_comm;
    parallel::mpi::Comm const                            distributed_particle_comm;

    static constexpr auto   tag        = parallel::mpi::Tag{ 875 };
    static constexpr char   null_dev[] = "/dev/null";
    static constexpr rank_t master{ 0 };
    [[nodiscard]] bool      is_subdomain_master() const { return master == subdomain_comm->rank(); }
    [[nodiscard]] bool      is_distributed_particle_master() const { return master == distributed_particle_comm.rank(); }

    // hdf5 space calculator
    template <class T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
    [[nodiscard]] static auto get_space(std::vector<T> const &payload)
    {
        auto mspace = hdf5::Space::simple(payload.size());
        mspace.select_all();
        auto fspace = hdf5::Space::simple(payload.size());
        fspace.select_all();
        return std::make_pair(std::move(mspace), std::move(fspace));
    }
    template <class T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
    [[nodiscard]] static auto get_space(std::vector<std::pair<T, T>> const &payload)
    {
        auto mspace = hdf5::Space::simple({ payload.size(), 2U });
        mspace.select_all();
        auto fspace = hdf5::Space::simple({ payload.size(), 2U });
        fspace.select_all();
        return std::make_pair(std::move(mspace), std::move(fspace));
    }
    [[nodiscard]] static auto get_space(std::vector<Scalar> const &payload) -> std::pair<hdf5::Space, hdf5::Space>;
    [[nodiscard]] static auto get_space(std::vector<Vector> const &payload) -> std::pair<hdf5::Space, hdf5::Space>;
    // exclude the off-diag components
    [[nodiscard]] static auto get_space(std::vector<Tensor> const &payload) -> std::pair<hdf5::Space, hdf5::Space>;
    [[nodiscard]] static auto get_space(std::vector<CurviCoord> const &payload) -> std::pair<hdf5::Space, hdf5::Space>;
    [[nodiscard]] static auto get_space(std::vector<Particle::PSD> const &payload) -> std::pair<hdf5::Space, hdf5::Space>;
};
PIC1D_END_NAMESPACE
