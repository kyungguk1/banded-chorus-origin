/*
 * Copyright (c) 2019-2021, Kyungguk Min
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

#include "BField.h"

#include "./EField.h"

P1D::BField::BField(ParamSet const &params) : GridQ{}, params{params}, geomtr{params}
{
    this->fill(geomtr.B0); // fill with background B
}

void P1D::BField::update(EField const &efield, Real const dt) noexcept
{
    Real const cdtODx = dt * params.c / params.Dx;
    _update(*this, efield, cdtODx);
}
void P1D::BField::_update(BField &B, EField const &E, Real const cdtODx) noexcept
{
    for (long i = 0; i < E.size(); ++i) {
        B[i].x += 0;
        B[i].y += (+E[i - 0].z - E[i - 1].z) * cdtODx;
        B[i].z += (-E[i - 0].y + E[i - 1].y) * cdtODx;
    }
}

namespace {
template <class Object>
decltype(auto) write_attr(Object &obj, [[maybe_unused]] P1D::BField const &bfield)
{
    return obj;
}
} // namespace
auto P1D::operator<<(hdf5::Group &obj, [[maybe_unused]] P1D::BField const &bfield) -> decltype(obj)
{
    return write_attr(obj, bfield);
}
auto P1D::operator<<(hdf5::Dataset &obj, [[maybe_unused]] P1D::BField const &bfield)
    -> decltype(obj)
{
    return write_attr(obj, bfield);
}
