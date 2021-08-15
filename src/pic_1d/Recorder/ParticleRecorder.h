/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include "Recorder.h"

#include <HDF5Kit/HDF5Kit.h>
#include <random>
#include <string>
#include <type_traits>
#include <vector>

PIC1D_BEGIN_NAMESPACE
/// marker particle recorder
/// field-aligned components are recorded;
/// suffix 1, 2, and 3 means three field-aligned components:
///     1 : parallel, 2 : perpendicular, and 3 : out-of-plane
///
class ParticleRecorder : public Recorder {
    std::mt19937 urbg;

public:
    explicit ParticleRecorder(parallel::mpi::Comm comm);

private:
    [[nodiscard]] std::string filepath(std::string const &wd, long step_count) const;

    void record(Domain const &domain, long step_count) override;
    void record_master(Domain const &domain, long step_count);
    void record_worker(Domain const &domain, long step_count);

    template <class Object>
    static decltype(auto) write_attr(Object &&obj, Domain const &domain, long step);
    static void           write_data(std::vector<Particle> ptls, hdf5::Group &root);
    template <class T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
    [[nodiscard]] static auto get_space(std::vector<T> const &payload);
    [[nodiscard]] static auto get_space(std::vector<Vector> const &payload);
    [[nodiscard]] static auto get_space(std::vector<Particle::PSD> const &payload);

    [[nodiscard]] std::vector<Particle> sample(PartSpecies const &sp, unsigned long max_count);
};
PIC1D_END_NAMESPACE
