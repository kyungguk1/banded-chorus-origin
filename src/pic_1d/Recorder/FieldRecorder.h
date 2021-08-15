/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include "Recorder.h"

#include <HDF5Kit/HDF5Kit.h>
#include <string>
#include <vector>

PIC1D_BEGIN_NAMESPACE
/// fluctuating (w/o background) electric and magnetic field recorder
/// field-aligned components are recorded;
/// suffix 1, 2, and 3 means three field-aligned components:
///     1 : parallel, 2 : perpendicular, and 3 : out-of-plane
///
class FieldRecorder : public Recorder {
public:
    explicit FieldRecorder(parallel::mpi::Comm comm);

private:
    [[nodiscard]] std::string filepath(std::string const &wd, long step_count) const;

    void record(Domain const &domain, long step_count) override;
    void record_master(Domain const &domain, long step_count);
    void record_worker(Domain const &domain, long step_count);

    template <class Object>
    static decltype(auto) write_attr(Object &&obj, Domain const &domain, long const step);
    static auto write_data(std::vector<Vector> payload, hdf5::Group &root, char const *name);

    [[nodiscard]] static std::vector<Vector> cart2fac(BField const &bfield, Geometry const &geomtr);
    [[nodiscard]] static std::vector<Vector> cart2fac(EField const &efield, Geometry const &geomtr);
};
PIC1D_END_NAMESPACE
