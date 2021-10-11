/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include "Recorder.h"

#include <fstream>
#include <string>

HYBRID1D_BEGIN_NAMESPACE
/// spatial average of field and ion energy density recorder
/// field-aligned components are recorded;
/// suffix 1, 2, and 3 means three field-aligned components:
///     1 : parallel, 2 : perpendicular, and 3 : out-of-plane
///
class EnergyRecorder : public Recorder {
    std::ofstream os;

public:
    EnergyRecorder(parallel::mpi::Comm subdomain_comm, parallel::mpi::Comm const &world_comm, ParamSet const &params);

private:
    [[nodiscard]] std::string filepath(std::string const &wd) const;

    void record(Domain const &domain, long step_count) override;

    [[nodiscard]] static Vector dump(BField const &bfield) noexcept;
    [[nodiscard]] static Vector dump(EField const &efield) noexcept;
    [[nodiscard]] static Tensor dump(Species const &sp) noexcept;
};
HYBRID1D_END_NAMESPACE
