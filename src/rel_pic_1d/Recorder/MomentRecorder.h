/*
 * Copyright (c) 2019-2022, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include "Recorder.h"

#include <string_view>

PIC1D_BEGIN_NAMESPACE
/// ion moment recorder
/// field-aligned components are recorded;
/// suffix 1, 2, and 3 means three field-aligned components:
///     1 : parallel, 2 : perpendicular, and 3 : out-of-plane
///
class MomentRecorder : public Recorder {
    [[nodiscard]] auto filepath(std::string_view const &wd, long step_count) const;

public:
    MomentRecorder(ParamSet const &params, parallel::mpi::Comm subdomain_comm, parallel::mpi::Comm const &world_comm);

private:
    void record(Domain const &domain, long step_count) override;
    void record_master(Domain const &domain, long step_count);
    void record_worker(Domain const &domain, long step_count);

    template <class Object>
    static decltype(auto) write_attr(Object &&obj, Domain const &domain, long step);
    template <class T>
    static auto write_data(std::vector<T> payload, hdf5::Group &root, char const *name);

public:
    [[nodiscard]] static auto cart_to_mfa(Grid<CartVector> const &mom1, Species const &) -> std::vector<MFAVector>;
    [[nodiscard]] static auto cart_to_mfa(Grid<FourCartTensor> const &mom2, Species const &) -> std::vector<FourMFATensor>;
};
PIC1D_END_NAMESPACE
