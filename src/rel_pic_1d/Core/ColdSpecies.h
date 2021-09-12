/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include "Species.h"

PIC1D_BEGIN_NAMESPACE
class EField;
class BField;

/// linearized cold fluid
///
class ColdSpecies : public Species {
    ColdPlasmaDesc desc;

public:
    [[nodiscard]] ColdPlasmaDesc const *operator->() const noexcept override { return &desc; }

    ScalarGrid mom0_full{}; // 0th moment on full grid
    VectorGrid mom1_full{}; // 1st moment on full grid

    ColdSpecies &operator=(ColdSpecies &&) = delete; // this should not be default-ed
    ColdSpecies(ParamSet const &params, ColdPlasmaDesc const &desc);
    ColdSpecies() = default; // needed for empty std::array

    void populate(); // load cold species; should only be called by master thread

    // update flow velocity by dt; <v>^n-1/2 -> <v>^n+1/2
    void update_vel(BField const &bfield, EField const &efield, Real dt);

    void collect_part(); // collect 0th & 1st moments
    void collect_all();  // collect all moments

private:
    void impl_update_nV(VectorGrid &nV, ScalarGrid const &n, Vector const &B0, EField const &E,
                        BorisPush const &boris) const;

    void impl_collect_part(ScalarGrid &n, VectorGrid &nV) const;
    void impl_collect_nuv(FourTensorGrid &nuv, ScalarGrid const &n, VectorGrid const &nV) const;
};
PIC1D_END_NAMESPACE
