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
#include <PIC/VT/Tensor.h>
#include <PIC/VT/Vector.h>

#include <cmath>
#include <type_traits>

LIBPIC_NAMESPACE_BEGIN(1)
namespace Detail {
template <class MirrorGeometry>
class CurviBasis {
    [[nodiscard]] inline decltype(auto) self() const noexcept { return static_cast<MirrorGeometry const &>(*this); }

    [[nodiscard]] inline static auto pow2(Real const x) noexcept { return x * x; }

    [[nodiscard]] CartVector covar_basis(CartCoord const &pos, std::integral_constant<long, 1>) const noexcept { return { self().D1() * (1 + self().xi2() * pow2(pos.x)), 0, 0 }; }
    [[nodiscard]] CartVector covar_basis(CartCoord const &pos, std::integral_constant<long, 2>) const noexcept { return { 0, self().D2() / std::sqrt(1 + self().xi2() * pow2(pos.x)), 0 }; }
    [[nodiscard]] CartVector covar_basis(CartCoord const &pos, std::integral_constant<long, 3>) const noexcept { return { 0, 0, self().D3() / std::sqrt(1 + self().xi2() * pow2(pos.x)) }; }
    [[nodiscard]] CartTensor covar_basis(CartCoord const &pos, std::integral_constant<long, 0>) const noexcept
    {
        auto const tmp      = 1 + self().xi2() * pow2(pos.x);
        auto const sqrt_tmp = std::sqrt(tmp);
        return { self().D1() * tmp, self().D2() / sqrt_tmp, self().D3() / sqrt_tmp, 0, 0, 0 };
    }

    [[nodiscard]] CartVector covar_basis(CurviCoord const &pos, std::integral_constant<long, 1>) const noexcept { return { self().D1() / pow2(std::cos(self().xi() * self().D1() * pos.q1)), 0, 0 }; }
    [[nodiscard]] CartVector covar_basis(CurviCoord const &pos, std::integral_constant<long, 2>) const noexcept { return { 0, self().D2() * std::cos(self().xi() * self().D1() * pos.q1), 0 }; }
    [[nodiscard]] CartVector covar_basis(CurviCoord const &pos, std::integral_constant<long, 3>) const noexcept { return { 0, 0, self().D3() * std::cos(self().xi() * self().D1() * pos.q1) }; }
    [[nodiscard]] CartTensor covar_basis(CurviCoord const &pos, std::integral_constant<long, 0>) const noexcept
    {
        auto const tmp = std::cos(self().xi() * self().D1() * pos.q1);
        return { self().D1() / pow2(tmp), self().D2() * tmp, self().D3() * tmp, 0, 0, 0 };
    }

    [[nodiscard]] CartVector contr_basis(CartCoord const &pos, std::integral_constant<long, 1>) const noexcept { return { self().inv_D1() / (1 + self().xi2() * pow2(pos.x)), 0, 0 }; }
    [[nodiscard]] CartVector contr_basis(CartCoord const &pos, std::integral_constant<long, 2>) const noexcept { return { 0, std::sqrt(1 + self().xi2() * pow2(pos.x)) * self().inv_D2(), 0 }; }
    [[nodiscard]] CartVector contr_basis(CartCoord const &pos, std::integral_constant<long, 3>) const noexcept { return { 0, 0, std::sqrt(1 + self().xi2() * pow2(pos.x)) * self().inv_D3() }; }
    [[nodiscard]] CartTensor contr_basis(CartCoord const &pos, std::integral_constant<long, 0>) const noexcept
    {
        auto const tmp      = 1 + self().xi2() * pow2(pos.x);
        auto const sqrt_tmp = std::sqrt(tmp);
        return { self().inv_D1() / tmp, sqrt_tmp * self().inv_D2(), sqrt_tmp * self().inv_D3(), 0, 0, 0 };
    }

    [[nodiscard]] CartVector contr_basis(CurviCoord const &pos, std::integral_constant<long, 1>) const noexcept { return { self().inv_D1() * pow2(std::cos(self().xi() * self().D1() * pos.q1)), 0, 0 }; }
    [[nodiscard]] CartVector contr_basis(CurviCoord const &pos, std::integral_constant<long, 2>) const noexcept { return { 0, self().inv_D2() / std::cos(self().xi() * self().D1() * pos.q1), 0 }; }
    [[nodiscard]] CartVector contr_basis(CurviCoord const &pos, std::integral_constant<long, 3>) const noexcept { return { 0, 0, self().inv_D3() / std::cos(self().xi() * self().D1() * pos.q1) }; }
    [[nodiscard]] CartTensor contr_basis(CurviCoord const &pos, std::integral_constant<long, 0>) const noexcept
    {
        auto const tmp = std::cos(self().xi() * self().D1() * pos.q1);
        return { self().inv_D1() * pow2(tmp), self().inv_D2() / tmp, self().inv_D3() / tmp, 0, 0, 0 };
    }

protected:
    CurviBasis() noexcept = default;

public:
    // MARK:- Metric Components
    /// Calculate covariant components of the metric tensor at a location
    ///
    [[nodiscard]] CovarTensor covar_metric(CartCoord const &pos) const noexcept
    {
        auto const tmp = 1 + self().xi2() * pow2(pos.x);
        return { pow2(self().D1() * tmp), pow2(self().D2()) / tmp, pow2(self().D3()) / tmp, 0, 0, 0 };
    }
    [[nodiscard]] CovarTensor covar_metric(CurviCoord const &pos) const noexcept
    {
        auto const tmp = pow2(std::cos(self().xi() * self().D1() * pos.q1));
        return { pow2(self().D1() / tmp), pow2(self().D2()) * tmp, pow2(self().D3()) * tmp, 0, 0, 0 };
    }

    /// Calculate contravariant components of the metric tensor at a location
    ///
    [[nodiscard]] ContrTensor contr_metric(CartCoord const &pos) const noexcept
    {
        auto const tmp = 1 + self().xi2() * pow2(pos.x);
        return { pow2(self().inv_D1() / tmp), tmp * pow2(self().inv_D2()), tmp * pow2(self().inv_D3()), 0, 0, 0 };
    }
    [[nodiscard]] ContrTensor contr_metric(CurviCoord const &pos) const noexcept
    {
        auto const tmp = pow2(std::cos(self().xi() * self().D1() * pos.q1));
        return { pow2(self().inv_D1() * tmp), pow2(self().inv_D2()) / tmp, pow2(self().inv_D3()) / tmp, 0, 0, 0 };
    }

    // MARK:- Basis
    /// Calculate i'th covariant basis vectors at a location
    /// \note The index '0' returns all three basis vectors.
    ///
    template <long i>
    [[nodiscard]] decltype(auto) covar_basis(CartCoord const &pos) const noexcept
    {
        static_assert(i >= 0 && i <= 3, "invalid index range");
        return covar_basis(pos, std::integral_constant<long, i>{});
    }
    template <long i>
    [[nodiscard]] decltype(auto) covar_basis(CurviCoord const &pos) const noexcept
    {
        static_assert(i >= 0 && i <= 3, "invalid index range");
        return covar_basis(pos, std::integral_constant<long, i>{});
    }

    /// Calculate i'th contravariant basis vectors at a location
    /// \note The index '0' returns all three basis vectors.
    ///
    template <long i>
    [[nodiscard]] decltype(auto) contr_basis(CartCoord const &pos) const noexcept
    {
        static_assert(i >= 0 && i <= 3, "invalid index range");
        return contr_basis(pos, std::integral_constant<long, i>{});
    }
    template <long i>
    [[nodiscard]] decltype(auto) contr_basis(CurviCoord const &pos) const noexcept
    {
        static_assert(i >= 0 && i <= 3, "invalid index range");
        return contr_basis(pos, std::integral_constant<long, i>{});
    }

    // MARK:- Vector Transformation
    /// Vector transformation from contravariant to covariant components
    /// \param contr Contravariant components of a vector.
    /// \param covar_metric Covariant components of the metric tensor at a given location.
    /// \return Covariant components of the transformed vector.
    [[nodiscard]] static decltype(auto) contr_to_covar(ContrVector const &contr, CovarTensor const &covar_metric) noexcept { return dot(contr, covar_metric); }

    /// Vector transformation from covariant to contravariant components
    /// \param covar Covariant components of a vector.
    /// \param contr_metric Contravariant components of the metric tensor at a given location.
    /// \return Contravariant components of the transformed vector.
    [[nodiscard]] static decltype(auto) covar_to_contr(CovarVector const &covar, ContrTensor const &contr_metric) noexcept { return dot(covar, contr_metric); }

    /// Vector transformation from Cartesian to contravariant components
    /// \param cart Cartesian components of a vector.
    /// \param contr_bases Three contravariant basis vectors at a given location.
    /// \return Contravariant components of the transformed vector.
    [[nodiscard]] static decltype(auto) cart_to_contr(CartVector const &cart, CartTensor const &contr_bases) noexcept { return ContrVector{ dot(contr_bases, cart) }; }

    /// Vector transformation from contravariant to Cartesian components
    /// \param contr Contravariant components of a vector.
    /// \param covar_bases Three covariant basis vectors at a given location.
    /// \return Cartesian components of the transformed vector.
    [[nodiscard]] static decltype(auto) contr_to_cart(ContrVector const &contr, CartTensor const &covar_bases) noexcept { return dot(CartVector{ contr }, covar_bases); }

    /// Vector transformation from Cartesian to covaraint components
    /// \param cart Cartesian components of a vector.
    /// \param covar_bases Three covariant basis vectors at a given location.
    /// \return Covariant components of the transformed vector.
    [[nodiscard]] static decltype(auto) cart_to_covar(CartVector const &cart, CartTensor const &covar_bases) noexcept { return CovarVector{ dot(covar_bases, cart) }; }

    /// Vector transformation from covariant to Cartesian components
    /// \param covar Covariant components of a vector.
    /// \param contr_bases Three contravariant basis vectors at a given location.
    /// \return Cartesian components of the transformed vector.
    [[nodiscard]] static decltype(auto) covar_to_cart(CovarVector const &covar, CartTensor const &contr_bases) noexcept { return dot(CartVector{ covar }, contr_bases); }

    /// Vector transformation from contravariant to covariant components
    /// \tparam Coord Coordinate type.
    /// \param contr Contravariant components of a vector.
    /// \param pos Current location.
    /// \return Covariant components of the transformed vector.
    template <class Coord>
    [[nodiscard]] decltype(auto) contr_to_covar(ContrVector const &contr, Coord const &pos) const noexcept { return contr_to_covar(contr, covar_metric(pos)); }

    /// Vector transformation from covariant to contravariant components
    /// \tparam Coord Coordinate type.
    /// \param covar Covariant components of a vector.
    /// \param pos Current location.
    /// \return Contravariant components of the transformed vector.
    template <class Coord>
    [[nodiscard]] decltype(auto) covar_to_contr(CovarVector const &covar, Coord const &pos) const noexcept { return covar_to_contr(covar, contr_metric(pos)); }

    /// Vector transformation from Cartesian to contravariant components
    /// \tparam Coord Coordinate type.
    /// \param cart Cartesian components of a vector.
    /// \param pos Current location.
    /// \return Contravariant components of the transformed vector.
    template <class Coord>
    [[nodiscard]] decltype(auto) cart_to_contr(CartVector const &cart, Coord const &pos) const noexcept { return cart_to_contr(cart, contr_basis<0>(pos)); }

    /// Vector transformation from contravariant to Cartesian components
    /// \tparam Coord Coordinate type.
    /// \param contr Contravariant components of a vector.
    /// \param pos Current location.
    /// \return Cartesian components of the transformed vector.
    template <class Coord>
    [[nodiscard]] decltype(auto) contr_to_cart(ContrVector const &contr, Coord const &pos) const noexcept { return contr_to_cart(contr, covar_basis<0>(pos)); }

    /// Vector transformation from Cartesian to covaraint components
    /// \tparam Coord Coordinate type.
    /// \param cart Cartesian components of a vector.
    /// \param pos Current location.
    /// \return Covariant components of the transformed vector.
    template <class Coord>
    [[nodiscard]] decltype(auto) cart_to_covar(CartVector const &cart, Coord const &pos) const noexcept { return cart_to_covar(cart, covar_basis<0>(pos)); }

    /// Vector transformation from covariant to Cartesian components
    /// \tparam Coord Coordinate type.
    /// \param covar Covariant components of a vector.
    /// \param pos Current location.
    /// \return Cartesian components of the transformed vector.
    template <class Coord>
    [[nodiscard]] decltype(auto) covar_to_cart(CovarVector const &covar, Coord const &pos) const noexcept { return covar_to_cart(covar, contr_basis<0>(pos)); }
};
} // namespace Detail
LIBPIC_NAMESPACE_END(1)
