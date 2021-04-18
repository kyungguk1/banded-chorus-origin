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

#include "MaxwellianVDF.h"

#include <cmath>

auto H1D::VDF::make(const BiMaxPlasmaDesc &desc) -> std::unique_ptr<VDF>
{
    return std::make_unique<MaxwellianVDF>(desc);
}

H1D::MaxwellianVDF::MaxwellianVDF(BiMaxPlasmaDesc const &desc) : VDF{}
{ // parameter check is assumed to be done already
    vth1  = std::sqrt(desc.beta1) * Input::c * std::abs(desc.Oc) / desc.op;
    T2OT1 = desc.T2_T1;
    xd    = desc.Vd / vth1;
    //
    vth1_cubed = vth1 * vth1 * vth1;
}

auto H1D::MaxwellianVDF::f0(Vector const &v) const noexcept -> Real
{
    // note that vel = {v1, v2, v3}/vth1
    // f0(x1, x2, x3) = exp(-(x1 - xd)^2)/√π * exp(-(x2^2 + x3^2)/(T2/T1))/(π T2/T1)
    //
    Real const x1_xd      = v.x - xd;
    Real const f1         = std::exp(-x1_xd * x1_xd) * M_2_SQRTPI * .5;
    Real const x2_squared = v.y * v.y + v.z * v.z;
    Real const f2         = std::exp(-x2_squared / T2OT1) / (M_PI * T2OT1);
    return f1 * f2;
}

auto H1D::MaxwellianVDF::variate() const -> Particle
{
    Particle ptl = load();

    // rescale
    //
    ptl.vel *= vth1;
    ptl.pos_x *= Input::Nx; // [0, Nx)

    // delta-f parameters
    //
    ptl.f = f0(ptl);
    // ptl.fOg = ptl.f/ptl.g0(ptl);
    static_assert(Particle::fOg == 1.0, "f and g should be identical");

    return ptl;
}
auto H1D::MaxwellianVDF::load() const -> Particle
{
    // position
    //
    Real const pos_x = bit_reversed<2>(); // [0, 1]

    // velocity in field-aligned frame (Hu et al., 2010, doi:10.1029/2009JA015158)
    //
    Real const phi1 = bit_reversed<3>() * 2 * M_PI;                               // [0, 2pi]
    Real const v1   = std::sqrt(-std::log(uniform_real<100>())) * std::sin(phi1); // v_para
    //
    Real const phi2 = bit_reversed<5>() * 2 * M_PI; // [0, 2pi]
    Real const _v2  = std::sqrt(-std::log(uniform_real<200>()) * T2OT1);
    Real const v2   = std::cos(phi2) * _v2; // in-plane v_perp
    Real const v3   = std::sin(phi2) * _v2; // out-of-plane v_perp

    // velocity in Cartesian frame
    //
    Vector const vel = geomtr.fac2cart({v1 + xd, v2, v3});

    return Particle{vel, pos_x};
}
