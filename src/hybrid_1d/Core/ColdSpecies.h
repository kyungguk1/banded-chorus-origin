/*
 * Copyright (c) 2019, Kyungguk Min
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * - Redistributions of source code must retain the above copyright notice, this
 *  list of conditions and the following disclaimer.
 *
 * - Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef ColdSpecies_h
#define ColdSpecies_h

#include "./Species.h"

HYBRID1D_BEGIN_NAMESPACE
class EField;
class BField;

/// linearized cold fluid
///
class ColdSpecies : public Species {
    ColdPlasmaDesc desc;

public:
    ScalarGrid mom0_full{}; // 0th moment on full grid
    VectorGrid mom1_full{}; // 1st moment on full grid
private:
    VectorGrid vect_buff{}; // vector buffer
public:
    [[nodiscard]] ColdPlasmaDesc const *operator->() const noexcept override { return &desc; }

    ColdSpecies &operator=(ColdSpecies const &) = default;
    ColdSpecies &operator=(ColdSpecies &&) = delete;

    ColdSpecies() = default; // needed for empty std::array
    explicit ColdSpecies(ParamSet const &params) : Species{ params } {} // needed for Domain_PC
    ColdSpecies(ParamSet const &params, ColdPlasmaDesc const &desc);

    void populate(); // load cold species; should only be called by master thread

    void update_den(Real dt); // update fluid number density by dt; <1>^n -> <1>^n+1
    void update_vel(BField const &bfield, EField const &efield,
                    Real dt); // update flow velocity by dt; <v>^n-1/2 -> <v>^n+1/2

    void collect_part(); // collect 0th & 1st moments
    void collect_all();  // collect all moments

private:
    void _update_n(ScalarGrid &n, VectorGrid const &nV, Real dt) const;
    void _update_nV(VectorGrid &new_nV, VectorGrid &old_nV, BorisPush pusher, ScalarGrid const &n,
                    VectorGrid const &B, EField const &E) const;

    void        _collect_part(ScalarGrid &n, VectorGrid &nV) const;
    static void _collect_nvv(TensorGrid &nvv, ScalarGrid const &n, VectorGrid const &nV);
};
HYBRID1D_END_NAMESPACE

#endif /* ColdSpecies_h */
