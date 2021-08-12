/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <Utility/Tensor.h>
#include <Utility/Vector.h>
#include <common-config.h>

#include <cmath>

COMMON_BEGIN_NAMESPACE
struct Geometry {
    using Real = double;

    // field-aligned unit vectors satisfying e1 = e2 x e3
    //
    Vector                  B0;               //!< the background magnetic field.
    Vector                  e1;               //!< parallel unit vector.
    Vector                  e2;               //!< in-plane perpendicular unit vector.
    static constexpr Vector e3 = { 0, 0, 1 }; //!< out-of-plane unit vector.

    Geometry() noexcept = default;

    /// Construct with the magnetic field vector
    ///
    /// \param O0 Background magnetic field represented as a reference cyclotron frequency.
    ///
    explicit Geometry(Vector const &O0);

    /// Construct with the magnetic field magnitude and an angle of the vector from the x-axis
    /// \details The magnetic field vector is assumed to lie in the x-y plane.
    ///
    /// \param O0 Magnitude of the background magnetic field represented as a reference cyclotron
    /// frequency.
    /// \param theta Angle between the magnetic field vector and the x-axis given in *radians*.
    ///
    Geometry(Real O0, Real theta) : Geometry({ O0 * std::cos(theta), O0 * std::sin(theta), 0 }) {}

    /// Transformation from field-aligned to cartesian
    ///
    /// \param v A field-aligned vector, {v1, v2, v3}.
    /// \return A cartesian vector, {vx, vy, vz}.
    ///
    [[nodiscard]] Vector fac2cart(Vector const &v) const noexcept
    {
        return e1 * v.x + e2 * v.y + e3 * v.z;
    }
    /// Transformation from field-aligned to cartesian
    ///
    /// \param vv A diagonal field-aligned tensor, {v1v1, v2v2, v3v3, 0, 0, 0}.
    /// \return A symmetric cartesian tensor, {vxvx, vyvy, vzvz, vxvy, vyvz, vzvx}.
    ///
    [[nodiscard]] Tensor fac2cart(Tensor const &vv) const noexcept
    {
        Tensor ret;
        ret.lo() = vv.xx * e1 * e1 + vv.yy * e2 * e2 + vv.zz * e3 * e3;
        ret.xy   = e1.x * e1.y * vv.xx + e2.x * e2.y * vv.yy + e3.x * e3.y * vv.zz;
        ret.yz   = e1.y * e1.z * vv.xx + e2.y * e2.z * vv.yy + e3.y * e3.z * vv.zz;
        ret.zx   = e1.x * e1.z * vv.xx + e2.x * e2.z * vv.yy + e3.x * e3.z * vv.zz;
        return ret;
    }

    /// Transformation from cartesian to field-aligned
    ///
    /// \param v A cartesian vector, {vx, vy, vz}.
    /// \return A field-aligned vector, {v1, v2, v3}.
    ///
    [[nodiscard]] Vector cart2fac(Vector const &v) const noexcept
    {
        return { dot(e1, v), dot(e2, v), dot(e3, v) };
    }
    /// Transformation from cartesian to field-aligned
    ///
    /// \param vv A symmetric cartesian tensor, {vxvx, vyvy, vzvz, vxvy, vyvz, vzvx}.
    /// \return A diagonal field-aligned tensor, {v1v1, v2v2, v3v3, 0, 0, 0}.
    ///
    [[nodiscard]] Tensor cart2fac(Tensor const &vv) const noexcept
    {
        return {
            dot(vv.lo(), e1 * e1)
                + 2 * (vv.xy * e1.x * e1.y + vv.zx * e1.x * e1.z + vv.yz * e1.y * e1.z),
            dot(vv.lo(), e2 * e2)
                + 2 * (vv.xy * e2.x * e2.y + vv.zx * e2.x * e2.z + vv.yz * e2.y * e2.z),
            dot(vv.lo(), e3 * e3)
                + 2 * (vv.xy * e3.x * e3.y + vv.zx * e3.x * e3.z + vv.yz * e3.y * e3.z),
            0,
            0,
            0,
        };
    }
};
COMMON_END_NAMESPACE
