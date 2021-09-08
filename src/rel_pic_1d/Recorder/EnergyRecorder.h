/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include "Recorder.h"

#include <fstream>
#include <string>

PIC1D_BEGIN_NAMESPACE
/// spatial average of field and ion energy density recorder
/// field-aligned components are recorded;
/// suffix 1, 2, and 3 means three field-aligned components:
///     1 : parallel, 2 : perpendicular, and 3 : out-of-plane
///
class EnergyRecorder : public Recorder {
    std::ofstream os;

public:
    EnergyRecorder(parallel::mpi::Comm comm, ParamSet const &params);

private:
    [[nodiscard]] std::string filepath(std::string const &wd) const;

    void record(Domain const &domain, long step_count) override;

    [[nodiscard]] static auto dump(BField const &bfield) noexcept -> Vector;
    [[nodiscard]] static auto dump(EField const &efield) noexcept -> Vector;
    [[nodiscard]] static auto dump(Species const &sp) noexcept -> FourTensor;
};
PIC1D_END_NAMESPACE
