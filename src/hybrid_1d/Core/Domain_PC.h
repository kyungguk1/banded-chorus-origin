/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include "Domain.h"

HYBRID1D_BEGIN_NAMESPACE
/// predictor-corrector algorithm
///
class Domain_PC : public Domain {
    BField      bfield_1;
    EField      efield_1;
    PartSpecies part_predict;
    ColdSpecies cold_predict;

public:
    explicit Domain_PC(ParamSet const &params, Delegate *delegate);

private:
    void advance_by(unsigned n_steps) override;
    void cycle(Domain const &domain);
    void predictor_step(Domain const &domain);
    void corrector_step(Domain const &domain);
};
HYBRID1D_END_NAMESPACE
