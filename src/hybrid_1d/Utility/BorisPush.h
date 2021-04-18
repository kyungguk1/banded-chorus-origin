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

#ifndef BorisPush_h
#define BorisPush_h

#include "../Macros.h"
#include "../Predefined.h"
#include "./Vector.h"

HYBRID1D_BEGIN_NAMESPACE
class BorisPush {
public:
    Real dt_2{};
    Real dtOc_2O0{};
    Real cDtOc_2O0{};

public:
    constexpr explicit BorisPush() noexcept = delete;
    constexpr explicit BorisPush(Real const dt, Real const c, Real const O0, Real const Oc) noexcept
    {
        dt_2      = 0.5 * dt;
        dtOc_2O0  = Oc * dt_2 / O0;
        cDtOc_2O0 = c * dtOc_2O0;
    }

    void operator()(Vector &V, Vector B, Vector cE, Real nu) const noexcept
    {
        nu *= dt_2;
        B *= dtOc_2O0;
        cE *= cDtOc_2O0;
        //
        // first half acceleration
        //
        V += (cE - nu * V) / (1 + nu / 2);
        //
        // rotation
        //
        V += rotate(V, B);
        //
        // second half acceleration
        //
        V += (cE - nu * V) / (1 + nu / 2);
    }
    void operator()(Vector &v, Vector B, Vector cE) const noexcept
    {
        B *= dtOc_2O0;
        cE *= cDtOc_2O0;
        //
        // first half acceleration
        //
        v += cE;
        //
        // rotation
        //
        v += rotate(v, B);
        //
        // second half acceleration
        //
        v += cE;
    }

private:
    [[nodiscard]] constexpr static Vector rotate(Vector const &v, Vector const &B) noexcept
    {
        return cross(v + cross(v, B), (2 / (1 + dot(B, B))) * B);
    }
};
HYBRID1D_END_NAMESPACE

#endif /* BorisPush_h */
