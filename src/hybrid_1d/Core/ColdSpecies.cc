/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "ColdSpecies.h"
#include "BField.h"
#include "EField.h"

#include <algorithm>

HYBRID1D_BEGIN_NAMESPACE
ColdSpecies::ColdSpecies(ParamSet const &params, ColdPlasmaDesc const &desc)
: Species{ params }, desc{ desc }
{
}
void ColdSpecies::populate()
{
    // initialize equilibrium moments
    //
    auto &n  = mom0_full;
    auto &nV = mom1_full;
    //
    constexpr Scalar n0{ 1 };
    Vector const     nV0 = Real{ n0 } * desc.Vd / params.O0 * params.geomtr.B0;
    for (long i = 0; i < nV.size(); ++i) { // only the interior
        n[i]  = n0;
        nV[i] = nV0;
    }
}

void ColdSpecies::update_vel(BField const &bfield, EField const &efield, Real const dt)
{
    impl_update_nV(mom1_full, mom0_full, bfield.params.geomtr.B0, efield,
                   BorisPush{ dt, params.c, params.O0, desc.Oc });
}
void ColdSpecies::impl_update_nV(VectorGrid &nV, ScalarGrid const &n, Vector const &B0,
                                 EField const &E, BorisPush const &boris) const
{
    for (long i = 0; i < nV.size(); ++i) {
        boris.non_relativistic(nV[i], B0, E[i] * Real{ n[i] });
    }
}

void ColdSpecies::collect_part()
{
    impl_collect_part(moment<0>(), moment<1>());
}
void ColdSpecies::collect_all()
{
    impl_collect_part(moment<0>(), moment<1>());
    impl_collect_nvv(moment<2>(), moment<0>(), moment<1>());
}
void ColdSpecies::impl_collect_part(ScalarGrid &n, VectorGrid &nV) const
{
    // must zero-out ghost cells
    //
    n.fill(Scalar{});
    std::copy(mom0_full.begin(), mom0_full.end(), n.begin());
    //
    nV.fill(Vector{});
    std::copy(mom1_full.begin(), mom1_full.end(), nV.begin());
}
void ColdSpecies::impl_collect_nvv(TensorGrid &nvv, ScalarGrid const &n, VectorGrid const &nV)
{
    for (long i = 0; i < nV.size(); ++i) {
        Tensor       &nvvi = nvv[i];
        Vector const &nVi  = nV[i];
        //
        nvvi.hi() = nvvi.lo()
            = nVi / Real{ n[i] };             // fill diag and off-diag terms with flow velocity
        nvvi.lo() *= nVi;                     // diagonal terms
        nvvi.hi() *= { nVi.y, nVi.z, nVi.x }; // off-diag terms
    }
}
HYBRID1D_END_NAMESPACE
