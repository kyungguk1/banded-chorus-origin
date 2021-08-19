/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include "Recorder.h"

#include <string>

PIC1D_BEGIN_NAMESPACE
/// ion moment recorder
/// field-aligned components are recorded;
/// suffix 1, 2, and 3 means three field-aligned components:
///     1 : parallel, 2 : perpendicular, and 3 : out-of-plane
///
class MomentRecorder : public Recorder {
public:
    explicit MomentRecorder(parallel::mpi::Comm comm);

private:
    [[nodiscard]] std::string filepath(std::string const &wd, long step_count) const;

    void record(Domain const &domain, long step_count) override;
    void record_master(Domain const &domain, long step_count);
    void record_worker(Domain const &domain, long step_count);

    template <class Object>
    static decltype(auto) write_attr(Object &&obj, Domain const &domain, long step);
    template <class T>
    static auto write_data(std::vector<T> payload, hdf5::Group &root, char const *name);

public:
    [[nodiscard]] static std::vector<Vector> cart2fac(VectorGrid const &mom1,
                                                      Geometry const   &geomtr);
    [[nodiscard]] static std::vector<Vector> cart2fac(TensorGrid const &mom2,
                                                      Geometry const   &geomtr);
};
PIC1D_END_NAMESPACE
