/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC/CartCoord.h>
#include <PIC/Config.h>
#include <PIC/CurviCoord.h>
#include <PIC/Predefined.h>
#include <PIC/Vector.h>

#include <cmath>

LIBPIC_BEGIN_NAMESPACE
namespace Detail {
template <class MirrorGeometry>
class MirrorField {
    [[nodiscard]] inline decltype(auto) self() const noexcept { return static_cast<MirrorGeometry const &>(*this); }

    [[nodiscard]] inline static auto pow2(Real const x) noexcept { return x * x; }
    [[nodiscard]] inline static auto pow4(Real const x) noexcept { return pow2(x) * pow2(x); }

protected:
    MirrorField() noexcept = default;

public:
    /// Contravariant components of B/B0
    ///
    [[nodiscard]] Vector Bcontr_div_B0(CartCoord const &) const noexcept { return { self().inv_D1(), 0, 0 }; }
    [[nodiscard]] Vector Bcontr_div_B0(CurviCoord const &) const noexcept { return { self().inv_D1(), 0, 0 }; }

    /// Covariant components of B/B0
    ///
    [[nodiscard]] Vector Bcovar_div_B0(CartCoord const &pos) const noexcept
    {
        return { self().D1() * pow2(1 + self().xi2() * pow2(pos.x)), 0, 0 };
    }
    [[nodiscard]] Vector Bcovar_div_B0(CurviCoord const &pos) const noexcept
    {
        return { self().D1() / pow4(std::cos(self().xi() * self().D1() * pos.q1)), 0, 0 };
    }

    /// Cartesian components of B/B0
    ///
    [[nodiscard]] Vector Bcart_div_B0(CartCoord const &pos) const noexcept
    {
        return { (1 + self().xi2() * pow2(pos.x)), 0, 0 };
    }
    /// Cartesian components of B/B0
    /// \param pos Cartesian x-component of position.
    /// \param pos_y Cartesian y-component of position.
    /// \param pos_z Cartesian z-component of position.
    ///
    [[nodiscard]] Vector Bcart_div_B0(CartCoord const &pos, Real pos_y, Real pos_z) const noexcept
    {
        auto const xi2 = self().xi2();
        return { (1 + xi2 * pos.x * pos.x), (0 - xi2 * pos.x * pos_y), (0 - xi2 * pos.x * pos_z) };
    }

    /// Cartesian components of B/B0
    ///
    [[nodiscard]] Vector Bcart_div_B0(CurviCoord const &pos) const noexcept
    {
        return { 1 / pow2(std::cos(self().xi() * self().D1() * pos.q1)), 0, 0 };
    }
    /// Cartesian components of B/B0
    /// \param pos Curvilinear q1-component of position.
    /// \param pos_y Cartesian y-component of position.
    /// \param pos_z Cartesian z-component of position.
    ///
    [[nodiscard]] Vector Bcart_div_B0(CurviCoord const &pos, Real pos_y, Real pos_z) const noexcept
    {
        return Bcart_div_B0(self().cotrans(pos), pos_y, pos_z);
    }

    /// Magnitude of B/B0
    ///
    [[nodiscard]] Real Bmag_div_B0(CartCoord const &pos) const noexcept { return 1 + self().xi2() * pow2(pos.x); }
    [[nodiscard]] Real Bmag_div_B0(CurviCoord const &pos) const noexcept { return 1 / pow2(std::cos(self().xi() * self().D1() * pos.q1)); }
};
} // namespace Detail
LIBPIC_END_NAMESPACE
