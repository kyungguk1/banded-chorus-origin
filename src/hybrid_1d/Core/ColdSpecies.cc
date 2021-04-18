//
// Copyright (c) 2019, Kyungguk Min
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// - Redistributions of source code must retain the above copyright notice, this
//  list of conditions and the following disclaimer.
//
// - Redistributions in binary form must reproduce the above copyright notice,
//  this list of conditions and the following disclaimer in the documentation
//  and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include "ColdSpecies.h"

#include "./BField.h"
#include "./EField.h"

#include <algorithm>
#include <stdexcept>

// helpers
//
namespace {
template <class T, long N>
auto const &full_grid(H1D::GridQ<T, N> &F, H1D::GridQ<T, N> const &H) noexcept
{
    for (long i = -H1D::Pad; i < F.size() + (H1D::Pad - 1); ++i) {
        (F[i] = H[i + 1] + H[i + 0]) *= 0.5;
    }
    return F;
}
//
constexpr bool enable_nonlinear_solver = false;
} // namespace

H1D::ColdSpecies::ColdSpecies(ParamSet const &params, ColdPlasmaDesc const &desc)
: Species{params}, desc{desc}
{
}
void H1D::ColdSpecies::populate()
{
    // initialize equilibrium moments
    //
    auto &n  = mom0_full;
    auto &nV = mom1_full;
    //
    constexpr Scalar n0{1};
    Vector const     nV0 = Real{n0} * desc.Vd / params.O0 * geomtr.B0;
    for (long i = 0; i < nV.size(); ++i) { // only the interior
        n[i]  = n0;
        nV[i] = nV0;
    }
}

void H1D::ColdSpecies::update_den(Real const dt)
{
    if constexpr (enable_nonlinear_solver)
        _update_n(mom0_full, mom1_full, dt);
}
void H1D::ColdSpecies::_update_n(ScalarGrid &n, VectorGrid const &nV, Real const dt) const
{
    static_assert(Pad >= 1, "not enough padding");
    for (long i = 0; i < n.size(); ++i) {
        Real const div_nV = (nV[i + 1].x - nV[i - 1].x) / (2 * params.Dx);
        n[i] -= dt * div_nV;
        if (Real{n[i]} < 0) {
            throw std::runtime_error{std::string{__FUNCTION__} + " - negative density"};
        }
    }
}

void H1D::ColdSpecies::update_vel(BField const &bfield, EField const &efield, Real const dt)
{
    moment<1>().fill(bfield.geomtr.B0);
    _update_nV(mom1_full, vect_buff = mom1_full, BorisPush{dt, params.c, params.O0, desc.Oc},
               mom0_full, moment<1>(), efield); // full_grid(moment<1>(), bfield), efield);
}
void H1D::ColdSpecies::_update_nV(VectorGrid &new_nV, VectorGrid &old_nV, BorisPush const pusher,
                                  ScalarGrid const &n, VectorGrid const &B, EField const &E) const
{
    static_assert(Pad >= 1, "not enough padding");

    // Lorentz force
    //
    for (long i = 0 - 1; i < new_nV.size() + 1; ++i) {
        pusher(new_nV[i], B[i], E[i] * Real{n[i]});
        (old_nV[i] += new_nV[i]) *= 0.5;
    }
    //
    // div nVV
    //
    for (long i = 0; enable_nonlinear_solver && i < new_nV.size(); ++i) {
        Vector const nVp1 = old_nV[i + 1], Vp1 = nVp1 / Real{n[i + 1]};
        Vector const nVm1 = old_nV[i - 1], Vm1 = nVm1 / Real{n[i - 1]};
        Vector const div_nVV = (nVp1.x * Vp1 - nVm1.x * Vm1) / (2 * params.Dx);
        new_nV[i] -= 2 * pusher.dt_2 * div_nVV;
    }
}

void H1D::ColdSpecies::collect_part()
{
    _collect_part(moment<0>(), moment<1>());
}
void H1D::ColdSpecies::collect_all()
{
    _collect_part(moment<0>(), moment<1>());
    _collect_nvv(moment<2>(), moment<0>(), moment<1>());
}
void H1D::ColdSpecies::_collect_part(ScalarGrid &n, VectorGrid &nV) const
{
    // must zero-out ghost cells
    //
    n.fill(Scalar{});
    std::copy(mom0_full.begin(), mom0_full.end(), n.begin());
    //
    nV.fill(Vector{});
    std::copy(mom1_full.begin(), mom1_full.end(), nV.begin());
}
void H1D::ColdSpecies::_collect_nvv(TensorGrid &nvv, ScalarGrid const &n, VectorGrid const &nV)
{
    for (long i = 0; i < nV.size(); ++i) {
        Tensor &      nvvi = nvv[i];
        Vector const &nVi  = nV[i];
        //
        nvvi.hi() = nvvi.lo() = nVi / Real{n[i]}; // fill diag and off-diag terms with flow velocity
        nvvi.lo() *= nVi;                         // diagonal terms
        nvvi.hi() *= {nVi.y, nVi.z, nVi.x};       // off-diag terms
    }
}