/*
 * Copyright (c) 2019, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#ifndef ColdSpecies_h
#define ColdSpecies_h

#include "./Species.h"

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

public:
    [[nodiscard]] ColdPlasmaDesc const *operator->() const noexcept override { return &desc; }

    ColdSpecies &operator=(ColdSpecies &&) = delete;
    ColdSpecies()                          = default; // needed for empty std::array
    ColdSpecies(ParamSet const &params, ColdPlasmaDesc const &desc);

    void populate(); // load cold species; should only be called by master thread

    void update_vel(BField const &bfield, EField const &efield,
                    Real dt); // update flow velocity by dt; <v>^n-1/2 -> <v>^n+1/2

    void collect_part(); // collect 0th & 1st moments
    void collect_all();  // collect all moments

private:
    void _update_nV(VectorGrid &nV, ScalarGrid const &n, Vector const &B0, EField const &E,
                    BorisPush pusher) const;

    void        _collect_part(ScalarGrid &n, VectorGrid &nV) const;
    static void _collect_nvv(TensorGrid &nvv, ScalarGrid const &n, VectorGrid const &nV);
};
PIC1D_END_NAMESPACE

#endif /* ColdSpecies_h */
