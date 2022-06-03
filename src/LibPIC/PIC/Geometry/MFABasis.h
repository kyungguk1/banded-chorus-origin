/*
 * Copyright (c) 2021-2022, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC/CartCoord.h>
#include <PIC/Config.h>
#include <PIC/CurviCoord.h>
#include <PIC/Predefined.h>
#include <PIC/VT/FourTensor.h>
#include <PIC/VT/FourVector.h>
#include <PIC/VT/Tensor.h>
#include <PIC/VT/Vector.h>

#include <cmath>

LIBPIC_NAMESPACE_BEGIN(1)
namespace Detail {
template <class MirrorGeometry>
class MFABasis {
    [[nodiscard]] inline decltype(auto) self() const noexcept { return static_cast<MirrorGeometry const &>(*this); }

    [[nodiscard]] inline static auto pow2(Real const x) noexcept { return x * x; }
    [[nodiscard]] inline static auto pow4(Real const x) noexcept { return pow2(x) * pow2(x); }

protected:
    MFABasis() noexcept = default;

public:
    /// Contravariant components of B/B0
    ///
    [[nodiscard]] ContrVector Bcontr_div_B0(CartCoord const &) const noexcept { return { self().inv_D1(), 0, 0 }; }
    [[nodiscard]] ContrVector Bcontr_div_B0(CurviCoord const &) const noexcept { return { self().inv_D1(), 0, 0 }; }

    /// Covariant components of B/B0
    ///
    [[nodiscard]] CovarVector Bcovar_div_B0(CartCoord const &pos) const noexcept
    {
        return { self().D1() * pow2(1 + self().xi2() * pow2(pos.x)), 0, 0 };
    }
    [[nodiscard]] CovarVector Bcovar_div_B0(CurviCoord const &pos) const noexcept
    {
        return { self().D1() / pow4(std::cos(self().xi() * self().D1() * pos.q1)), 0, 0 };
    }

    /// Cartesian components of B/B0
    ///
    [[nodiscard]] CartVector Bcart_div_B0(CartCoord const &pos) const noexcept
    {
        return { 1 + self().xi2() * pow2(pos.x), 0, 0 };
    }
    /// Cartesian components of B/B0
    /// \param pos Cartesian x-component of position.
    /// \param pos_y Cartesian y-component of position.
    /// \param pos_z Cartesian z-component of position.
    ///
    [[nodiscard]] CartVector Bcart_div_B0(CartCoord const &pos, Real pos_y, Real pos_z) const noexcept
    {
        auto const xi2 = self().xi2();
        return { 1 + xi2 * pos.x * pos.x, 0 - xi2 * pos.x * pos_y, 0 - xi2 * pos.x * pos_z };
    }

    /// Cartesian components of B/B0
    ///
    [[nodiscard]] CartVector Bcart_div_B0(CurviCoord const &pos) const noexcept
    {
        return { 1 / pow2(std::cos(self().xi() * self().D1() * pos.q1)), 0, 0 };
    }
    /// Cartesian components of B/B0
    /// \param pos Curvilinear q1-component of position.
    /// \param pos_y Cartesian y-component of position.
    /// \param pos_z Cartesian z-component of position.
    ///
    [[nodiscard]] CartVector Bcart_div_B0(CurviCoord const &pos, Real pos_y, Real pos_z) const noexcept
    {
        return Bcart_div_B0(self().cotrans(pos), pos_y, pos_z);
    }

    /// Magnitude of B/B0
    ///
    [[nodiscard]] Real Bmag_div_B0(CartCoord const &pos) const noexcept { return 1 + self().xi2() * pow2(pos.x); }
    [[nodiscard]] Real Bmag_div_B0(CurviCoord const &pos) const noexcept { return 1 / pow2(std::cos(self().xi() * self().D1() * pos.q1)); }

    // MARK:- Basis
    template <class Coord>
    [[nodiscard]] static decltype(auto) e1(Coord const &) noexcept { return CartVector{ 1, 0, 0 }; }
    template <class Coord>
    [[nodiscard]] static decltype(auto) e2(Coord const &) noexcept { return CartVector{ 0, 1, 0 }; }
    template <class Coord>
    [[nodiscard]] static decltype(auto) e3(Coord const &) noexcept { return CartVector{ 0, 0, 1 }; }

    // MARK:- Vector Transform
    // for the present 1D situation, all mfa <-> cart vector transformations at the central field line are just pass-through
    template <class Coord>
    [[nodiscard]] static decltype(auto) mfa_to_cart(MFAVector const &vmfa, Coord const &) noexcept { return CartVector{ vmfa }; }
    template <class Coord>
    [[nodiscard]] static decltype(auto) mfa_to_cart(MFATensor const &Tmfa, Coord const &) noexcept { return CartTensor{ Tmfa }; }

    template <class Coord>
    [[nodiscard]] static decltype(auto) cart_to_mfa(CartVector const &vcart, Coord const &) noexcept { return MFAVector{ vcart }; }
    template <class Coord>
    [[nodiscard]] static decltype(auto) cart_to_mfa(CartTensor const &Tcart, Coord const &) noexcept { return MFATensor{ Tcart }; }

    template <class Coord>
    [[nodiscard]] static decltype(auto) mfa_to_cart(FourMFAVector const &vmfa, Coord const &) noexcept { return FourCartVector{ vmfa }; }
    template <class Coord>
    [[nodiscard]] static decltype(auto) mfa_to_cart(FourMFATensor const &Tmfa, Coord const &) noexcept { return FourCartTensor{ Tmfa }; }

    template <class Coord>
    [[nodiscard]] static decltype(auto) cart_to_mfa(FourCartVector const &vcart, Coord const &) noexcept { return FourMFAVector{ vcart }; }
    template <class Coord>
    [[nodiscard]] static decltype(auto) cart_to_mfa(FourCartTensor const &vvcart, Coord const &) noexcept { return FourMFATensor{ vvcart }; }
};
} // namespace Detail
LIBPIC_NAMESPACE_END(1)
