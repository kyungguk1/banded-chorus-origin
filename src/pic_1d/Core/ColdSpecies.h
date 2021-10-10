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
    ScalarGrid mom0_full{}; // 0th moment on full grid
    VectorGrid mom1_full{}; // 1st moment on full grid

    [[nodiscard]] ColdPlasmaDesc const *operator->() const noexcept override { return &desc; }

    ColdSpecies &operator=(ColdSpecies &&) = delete; // this should not be default-ed
    ColdSpecies(ParamSet const &params, ColdPlasmaDesc const &desc);
    ColdSpecies() = default; // needed for empty std::array

    /// Load cold species
    /// \note This should only be called by master thread.
    /// \param color This is unused here; just to keep the symmetry with `PartSpecies::populate`.
    /// \param divisor The number of groups to which cold fluid are divided.
    void                populate(long color, long divisor);

    // update flow velocity by dt; <v>^n-1/2 -> <v>^n+1/2
    void update_vel(BField const &bfield, EField const &efield, Real dt);

    void collect_part(); // collect 0th & 1st moments
    void collect_all();  // collect all moments

private:
    void impl_update_nV(VectorGrid &nV, ScalarGrid const &n, EField const &E, BorisPush const &boris) const;

    void        impl_collect_part(ScalarGrid &n, VectorGrid &nV) const;
    static void impl_collect_nvv(TensorGrid &nvv, ScalarGrid const &n, VectorGrid const &nV);
};
PIC1D_END_NAMESPACE
