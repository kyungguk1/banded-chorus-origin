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
void ColdSpecies::populate(long, long const divisor)
{
    if (divisor <= 0)
        throw std::invalid_argument{ std::string{ __PRETTY_FUNCTION__ } + " - non-positive divisor" };

    // initialize equilibrium moments
    auto const   n0 = Scalar{ 1 } / divisor;
    Vector const V0 = { 0, 0, 0 };
    for (long i = 0; i < mom1_full.size(); ++i) { // only the interior
        mom0_full[i] = n0;
        mom1_full[i] = V0 * Real{ n0 };
    }
}

void ColdSpecies::update_vel(BField const &, EField const &efield, Real const dt)
{
    impl_update_nV(mom1_full, mom0_full, efield, BorisPush{ dt, params.c, params.O0, desc.Oc });
}
void ColdSpecies::impl_update_nV(VectorGrid &nV, ScalarGrid const &n, EField const &E, BorisPush const &boris) const
{
    auto const q1min = params.half_grid_subdomain_extent.min();
    for (long i = 0; i < nV.size(); ++i) {
        boris.non_relativistic(nV[i], geomtr.Bcart(CurviCoord{ i + q1min }), E[i] * Real{ n[i] });
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
