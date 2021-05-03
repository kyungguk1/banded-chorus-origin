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

#ifndef MaxwellianVDF_h
#define MaxwellianVDF_h

#include "./VDF.h"

HYBRID1D_BEGIN_NAMESPACE
/// bi-Maxwellian velocity distribution function
///
/// f(v1, v2) = exp(-(x1 - xd)^2 -x2^2)/(π^3/2 vth1^3 T2/T1),
/// where x1 = v1/vth1, xd = vd/vth1, x2 = v2/(vth1*√(T2/T1))), and
/// T2 and T1 are temperatures in directions perpendicular and
/// parallel to the background magnetic field direction, respectively
///
class MaxwellianVDF final : public VDF {
    Real vth1;  //!< Parallel thermal speed.
    Real T2OT1; //!< Temperature anisotropy, T2/T1.
    Real xd;    //!< Parallel drift speed normalized to vth1.
    //
    Real vth1_cubed;

public:
    explicit MaxwellianVDF(BiMaxPlasmaDesc const &desc);

public:
    [[nodiscard]] Scalar n0(Real) const override
    {
        constexpr Real n0 = 1;
        return n0;
    }
    [[nodiscard]] Vector nV0(Real const pos_x) const override
    {
        return geomtr.fac2cart({ xd * vth1, 0, 0 }) * Real{ n0(pos_x) };
    }
    [[nodiscard]] Tensor nvv0(Real const pos_x) const override
    {
        Tensor vv{ 1 + 2 * xd * xd, T2OT1, T2OT1, 0, 0, 0 }; // field-aligned 2nd moment
        return geomtr.fac2cart(vv *= .5 * vth1 * vth1) * Real{ n0(pos_x) };
    }

    [[nodiscard]] Real delta_f(Particle const &ptl) const override { return 1 - f0(ptl) / ptl.f; }

    [[nodiscard]] Particle variate() const override;

private:
    [[nodiscard]] Particle load() const;

    // equilibrium physical distribution function
    // f0(x1, x2, x3) = exp(-(x1 - xd)^2)/√π * exp(-(x2^2 + x3^2)/(T2/T1))/(π T2/T1)
    //
    [[nodiscard]] Real f0(Vector const &vel) const noexcept;
    [[nodiscard]] Real f0(Particle const &ptl) const noexcept
    {
        return f0(geomtr.cart2fac(ptl.vel) / vth1) / vth1_cubed;
    }

    // marker particle distribution function
    //
    [[nodiscard]] Real g0(Vector const &vel) const noexcept { return f0(vel); }
    [[nodiscard]] Real g0(Particle const &ptl) const noexcept
    {
        return g0(geomtr.cart2fac(ptl.vel) / vth1) / vth1_cubed;
    }
};
HYBRID1D_END_NAMESPACE

#endif /* MaxwellianVDF_h */
