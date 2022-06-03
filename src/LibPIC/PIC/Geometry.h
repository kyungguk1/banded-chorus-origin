/*
 * Copyright (c) 2019-2022, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC/Config.h>
#include <PIC/Geometry/MirrorGeometry.h>
#include <PIC/Predefined.h>
#include <PIC/VT/Vector.h>

#include <cmath>

LIBPIC_NAMESPACE_BEGIN(1)
class Geometry : public Detail::MirrorGeometry {
    Real m_O0{};

public:
    Geometry() noexcept = default;
    Geometry(Real xi, Real D1, Real O0)
    : MirrorGeometry{ xi, D1 }, m_O0{ O0 } {}

    /// Magnitude of B at the origin
    [[nodiscard]] auto B0() const noexcept { return m_O0; }

    /// Contravariant components of B
    template <class Coord>
    [[nodiscard]] decltype(auto) Bcontr(Coord const &pos) const noexcept { return Bcontr_div_B0(pos) * B0(); }

    /// Covariant components of B
    template <class Coord>
    [[nodiscard]] decltype(auto) Bcovar(Coord const &pos) const noexcept { return Bcovar_div_B0(pos) * B0(); }

    /// Cartesian components of B
    template <class Coord>
    [[nodiscard]] decltype(auto) Bcart(Coord const &pos) const noexcept { return Bcart_div_B0(pos) * B0(); }

    /// Cartesian components of B
    /// \tparam Coord Coordinate type.
    /// \param pos First component of position coordinates.
    /// \param pos_y Cartesian y-component of position.
    /// \param pos_z Cartesian z-component of position.
    template <class Coord>
    [[nodiscard]] decltype(auto) Bcart(Coord const &pos, Real pos_y, Real pos_z) const noexcept { return Bcart_div_B0(pos, pos_y, pos_z) * B0(); }

    /// Magnitude of B/B0
    template <class Coord>
    [[nodiscard]] Real Bmag(Coord const &pos) const noexcept { return Bmag_div_B0(pos) * B0(); }
};
LIBPIC_NAMESPACE_END(1)
