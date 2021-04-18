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

#include "EField.h"

#include "./BField.h"
#include "./Current.h"

P1D::EField::EField(ParamSet const &params) : GridQ{}, params{params}, geomtr{params}
{
}

void P1D::EField::update(BField const &bfield, Current const &current, Real const dt) noexcept
{
    Real const cdtODx = dt * params.c / params.Dx;
    _update(*this, bfield, cdtODx, current, dt);
}

void P1D::EField::_update(EField &E, BField const &B, Real const cdtODx, Current const &J,
                          Real const dt) noexcept
{
    for (long i = 0; i < B.size(); ++i) {
        E[i].x += 0;
        E[i].y += (-B[i + 1].z + B[i + 0].z) * cdtODx;
        E[i].z += (+B[i + 1].y - B[i + 0].y) * cdtODx;
        //
        E[i] -= J[i] * dt;
    }
}
