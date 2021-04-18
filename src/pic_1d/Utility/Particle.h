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

#ifndef Particle_h
#define Particle_h

#include "../Macros.h"
#include "../Predefined.h"
#include "./Vector.h"

#include <limits>
#include <ostream>

PIC1D_BEGIN_NAMESPACE
/// single particle description
///
struct Particle {
    static constexpr Real quiet_nan = std::numeric_limits<Real>::quiet_NaN();

    Vector vel{quiet_nan};   //!< 3-component velocity vector
    Real   pos_x{quiet_nan}; //!< x-component of position

    explicit Particle() noexcept = default;
    explicit Particle(Vector const &vel, Real const pos_x) noexcept : vel{vel}, pos_x{pos_x} {}

    // for delta-f
    //
    static constexpr Real fOg{1}; // f(0, x(0), v(0))/g(0, x(0), v(0)),
                                  // where g is the marker particle distribution
    Real f{quiet_nan};            // f(0, x(0), v(0))
    Real w{quiet_nan}; // f(0, x(0), v(0))/g(0, x(0), v(0)) - f_0(x(t), v(t))/g(0, x(0), v(0))

    // pretty print
    //
    template <class CharT, class Traits>
    friend decltype(auto) operator<<(std::basic_ostream<CharT, Traits> &os, Particle const &ptl)
    {
        return os << '{' << ptl.vel << ", " << '{' << ptl.pos_x << '}' << ", " << '{' << ptl.f
                  << ", " << ptl.w << '}' << '}';
    }
};
PIC1D_END_NAMESPACE

#endif /* Particle_h */
