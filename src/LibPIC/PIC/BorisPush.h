/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC/Config.h>
#include <PIC/Predefined.h>
#include <PIC/Vector.h>

LIBPIC_BEGIN_NAMESPACE
struct BorisPush {
    Real c2;        // c^2
    Real dt_2;      // dt/2
    Real dtOc_2O0;  // (dt/2) * (Oc/O0)
    Real cDtOc_2O0; // c * (dt/2) * (Oc/O0)

    BorisPush(Real const dt, Real const c, Real const O0, Real const Oc) noexcept
    {
        c2        = c * c;
        dt_2      = 0.5 * dt;
        dtOc_2O0  = Oc * dt_2 / O0;
        cDtOc_2O0 = c * dtOc_2O0;
    }

    [[deprecated]] void resistive(Vector &V, Vector B, Vector E, Real nu) const noexcept
    {
        nu *= dt_2;
        B *= dtOc_2O0;
        auto const &cE = E *= cDtOc_2O0;
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

    /// Non-relativistic Boris push
    ///
    /// @param [in,out] v Particle's velocity
    /// @param B Magnetic field at particle's position
    /// @param E Electric field at particle's position
    ///
    void non_relativistic(Vector &v, Vector B, Vector E) const noexcept
    {
        B *= dtOc_2O0;
        auto const &cE = E *= cDtOc_2O0;
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

    /// Relativistic Boris push
    ///
    /// @param [in,out] gv gamma * v, i.e., relativistic momentum
    /// @param B Magnetic field at particle's position
    /// @param E Electric field at particle's position
    /// @return Updated relativistic factor
    ///
    [[nodiscard]] Real relativistic(Vector &gv, Vector B, Vector E) const noexcept
    {
        B *= dtOc_2O0;
        auto const &cE = E *= cDtOc_2O0;

        // first half acceleration
        gv += cE;

        // rotation
        gv += rotate(gv, B /= gamma(gv));

        // second half acceleration
        gv += cE;

        return gamma(gv);
    }

private:
    [[nodiscard]] static Vector rotate(Vector const &v, Vector const &B) noexcept
    {
        return cross(v + cross(v, B), (2 / (1 + dot(B, B))) * B);
    }
    [[nodiscard]] Real gamma(Vector const &gv) const noexcept
    {
        return std::sqrt(1 + dot(gv, gv) / c2);
    }
};
LIBPIC_END_NAMESPACE
