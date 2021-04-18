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

#ifndef Geometry_h
#define Geometry_h

#include "./InputWrapper.h"
#include "./Utility/Vector.h"

HYBRID1D_BEGIN_NAMESPACE
class [[nodiscard]] Geometry {
public:
    // field-aligned unit vectors satisfying e1 = e2 x e3
    //
    static constexpr Vector e3 = {0, 0, 1}; //!< out-of-plane unit vector.
    Vector                  e1;             //!< parallel unit vector.
    Vector                  e2;             //!< in-plane perpendicular unit vector.
    Vector                  B0;             //!< the background magnetic field.

public:
    explicit Geometry(Input const &params) noexcept;

public:
    // transformation from field-aligned to cartesian
    //
    [[nodiscard]] Vector fac2cart(Vector const &v) const noexcept
    {
        // v = {v1, v2, v3}
        //
        return e1 * v.x + e2 * v.y + e3 * v.z;
    }
    [[nodiscard]] Tensor fac2cart(Tensor const &vv) const noexcept
    {
        // vv = {v1v1, v2v2, v3v3, 0, 0, 0}
        //
        Tensor ret;
        ret.lo() = vv.xx * e1 * e1 + vv.yy * e2 * e2 + vv.zz * e3 * e3;
        ret.xy   = e1.x * e1.y * vv.xx + e2.x * e2.y * vv.yy + e3.x * e3.y * vv.zz;
        ret.yz   = e1.y * e1.z * vv.xx + e2.y * e2.z * vv.yy + e3.y * e3.z * vv.zz;
        ret.zx   = e1.x * e1.z * vv.xx + e2.x * e2.z * vv.yy + e3.x * e3.z * vv.zz;
        return ret;
    }

    // transformation from cartesian to field-aligned
    //
    [[nodiscard]] Vector cart2fac(Vector const &v) const noexcept
    {
        // {v_||, v_perp, v_z}
        //
        return {dot(e1, v), dot(e2, v), dot(e3, v)};
    }
    [[nodiscard]] Vector cart2fac(Tensor const &vv) const noexcept
    {
        // similarity transformation
        // {v_||^2, v_perp^2, v_z^2}
        //
        return {dot(vv.lo(), e1 * e1)
                    + 2 * (vv.xy * e1.x * e1.y + vv.zx * e1.x * e1.z + vv.yz * e1.y * e1.z),
                dot(vv.lo(), e2 * e2)
                    + 2 * (vv.xy * e2.x * e2.y + vv.zx * e2.x * e2.z + vv.yz * e2.y * e2.z),
                dot(vv.lo(), e3 * e3)
                    + 2 * (vv.xy * e3.x * e3.y + vv.zx * e3.x * e3.z + vv.yz * e3.y * e3.z)};
    }
};
HYBRID1D_END_NAMESPACE

#endif /* Geometry_h */
