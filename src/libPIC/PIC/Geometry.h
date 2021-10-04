/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC/Config.h>
#include <PIC/FourTensor.h>
#include <PIC/FourVector.h>
#include <PIC/MirrorGeometry.h>
#include <PIC/Predefined.h>
#include <PIC/Tensor.h>
#include <PIC/Vector.h>

#include <cmath>

LIBPIC_BEGIN_NAMESPACE
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
    [[nodiscard]] Vector Bcontr(Coord const &pos) const noexcept { return Bcontr_div_B0(pos) * B0(); }

    /// Covariant components of B
    template <class Coord>
    [[nodiscard]] Vector Bcovar(Coord const &pos) const noexcept { return Bcovar_div_B0(pos) * B0(); }

    /// Cartesian components of B
    template <class Coord>
    [[nodiscard]] Vector Bcart(Coord const &pos) const noexcept { return Bcart_div_B0(pos) * B0(); }

    /// Cartesian components of B
    /// \tparam Coord Coordinate type.
    /// \param pos First component of position coordinates.
    /// \param pos_y Cartesian y-component of position.
    /// \param pos_z Cartesian z-component of position.
    template <class Coord>
    [[nodiscard]] Vector Bcart(Coord const &pos, Real pos_y, Real pos_z) const noexcept { return Bcart_div_B0(pos, pos_y, pos_z) * B0(); }

    /// Magnitude of B/B0
    template <class Coord>
    [[nodiscard]] Real Bmag(Coord const &pos) const noexcept { return Bmag_div_B0(pos) * B0(); }

    // for the present 1D situation, all fac <-> cart vector transformations at the central field line are just pass-through
    template <class Coord>
    [[nodiscard]] static decltype(auto) fac_to_cart(Vector const &vfac, Coord const &) noexcept { return vfac; }
    template <class Coord>
    [[nodiscard]] static decltype(auto) fac_to_cart(Tensor const &vvfac, Coord const &) noexcept { return vvfac; }

    template <class Coord>
    [[nodiscard]] static decltype(auto) cart_to_fac(Vector const &vcart, Coord const &) noexcept { return vcart; }
    template <class Coord>
    [[nodiscard]] static decltype(auto) cart_to_fac(Tensor const &vvcart, Coord const &) noexcept { return vvcart; }

    template <class Coord>
    [[nodiscard]] static decltype(auto) fac_to_cart(FourVector const &vfac, Coord const &) noexcept { return vfac; }
    template <class Coord>
    [[nodiscard]] static decltype(auto) fac_to_cart(FourTensor const &vvfac, Coord const &) noexcept { return vvfac; }

    template <class Coord>
    [[nodiscard]] static decltype(auto) cart_to_fac(FourVector const &vcart, Coord const &) noexcept { return vcart; }
    template <class Coord>
    [[nodiscard]] static decltype(auto) cart_to_fac(FourTensor const &vvcart, Coord const &) noexcept { return vvcart; }

    // field-aligned vectors at the central field line
    template <class Coord>
    [[nodiscard]] static constexpr Vector e1(Coord const &) noexcept { return { 1, 0, 0 }; }
    template <class Coord>
    [[nodiscard]] static constexpr Vector e2(Coord const &) noexcept { return { 0, 1, 0 }; }
    template <class Coord>
    [[nodiscard]] static constexpr Vector e3(Coord const &) noexcept { return { 0, 0, 1 }; }
};
LIBPIC_END_NAMESPACE
