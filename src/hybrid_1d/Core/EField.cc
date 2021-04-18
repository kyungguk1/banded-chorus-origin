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

#include "EField.h"

#include "./BField.h"
#include "./Charge.h"
#include "./Current.h"

#include <cmath>

H1D::EField::EField(ParamSet const &params) : GridQ{}, Je{}, Pe{}, params{params}, geomtr{params}
{
}

void H1D::EField::update(BField const &bfield, Charge const &charge,
                         Current const &current) noexcept
{
    _update_Pe(Pe, charge);
    _update_Je(Je, current, bfield);
    _update_E(*this, bfield, charge);
}

void H1D::EField::_update_Pe(ScalarGrid &Pe, Charge const &rho) const noexcept
{
    eFluidDesc const ef = params.efluid_desc;
    //
    Real const O0        = params.O0;
    Real const O02beO2   = (O0 * O0) * ef.beta * 0.5;
    Real const mOeOO0oe2 = -ef.Oc / (O0 * (ef.op * ef.op));
    for (long i = -Pad; i < Pe.size() + Pad; ++i) {
        Pe[i] = std::pow(mOeOO0oe2 * Real{rho[i]}, ef.gamma) * O02beO2;
    }
}
void H1D::EField::_update_Je(VectorGrid &Je, Current const &Ji, BField const &B) const noexcept
{
    Real const cODx = params.c / params.Dx;
    for (long i = 0; i < B.size(); ++i) {
        // J total
        //
        Je[i].x = 0;
        Je[i].y = (-B[i + 1].z + B[i + 0].z) * cODx;
        Je[i].z = (+B[i + 1].y - B[i + 0].y) * cODx;

        // Je = J - Ji
        //
        Je[i] -= Ji[i];
    }
}
void H1D::EField::_update_E(EField &E, BField const &B, Charge const &rho) const noexcept
{
    Real const cODx = params.c / params.Dx;
    for (long i = 0; i < E.size(); ++i) {
        Vector &Ei = E[i];

        // 1. Je x B term
        //
        Ei = cross(Je[i], (B[i + 1] + B[i + 0]) * 0.5);

        // 2. pressure gradient term
        //
        Ei.x -= 0.5 * Real{Pe[i + 1] - Pe[i - 1]} * cODx;
        Ei.y -= 0;
        Ei.z -= 0;

        // 3. divide by charge density
        //
        Ei /= Real{rho[i]};
    }
}
