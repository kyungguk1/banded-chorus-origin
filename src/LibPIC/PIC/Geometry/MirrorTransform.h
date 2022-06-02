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

LIBPIC_NAMESPACE_BEGIN(1)
namespace Detail {
template <class MirrorGeometry>
class MirrorTransform {
    [[nodiscard]] inline decltype(auto) self() const noexcept { return static_cast<MirrorGeometry const &>(*this); }

protected:
    MirrorTransform() noexcept = default;

public:
    /// Vector transformation from contravariant to covariant components
    /// \param contr Contravariant components of a vector.
    /// \param covar_metric Covariant components of the metric tensor at a given location.
    /// \return Covariant components of the transformed vector.
    [[nodiscard]] static Vector contr_to_covar(Vector const &contr, Tensor const &covar_metric) noexcept { return dot(contr, covar_metric); }

    /// Vector transformation from covariant to contravariant components
    /// \param covar Covariant components of a vector.
    /// \param contr_metric Contravariant components of the metric tensor at a given location.
    /// \return Contravariant components of the transformed vector.
    [[nodiscard]] static Vector covar_to_contr(Vector const &covar, Tensor const &contr_metric) noexcept { return dot(covar, contr_metric); }

    /// Vector transformation from Cartesian to contravariant components
    /// \param cart Cartesian components of a vector.
    /// \param contr_bases Three contravariant basis vectors at a given location.
    /// \return Contravariant components of the transformed vector.
    [[nodiscard]] static Vector cart_to_contr(Vector const &cart, Tensor const &contr_bases) noexcept { return dot(contr_bases, cart); }

    /// Vector transformation from contravariant to Cartesian components
    /// \param contr Contravariant components of a vector.
    /// \param covar_bases Three covariant basis vectors at a given location.
    /// \return Cartesian components of the transformed vector.
    [[nodiscard]] static Vector contr_to_cart(Vector const &contr, Tensor const &covar_bases) noexcept { return dot(contr, covar_bases); }

    /// Vector transformation from Cartesian to covaraint components
    /// \param cart Cartesian components of a vector.
    /// \param covar_bases Three covariant basis vectors at a given location.
    /// \return Covariant components of the transformed vector.
    [[nodiscard]] static Vector cart_to_covar(Vector const &cart, Tensor const &covar_bases) noexcept { return dot(covar_bases, cart); }

    /// Vector transformation from covariant to Cartesian components
    /// \param covar Covariant components of a vector.
    /// \param contr_bases Three contravariant basis vectors at a given location.
    /// \return Cartesian components of the transformed vector.
    [[nodiscard]] static Vector covar_to_cart(Vector const &covar, Tensor const &contr_bases) noexcept { return dot(covar, contr_bases); }

    /// Vector transformation from contravariant to covariant components
    /// \tparam Coord Coordinate type.
    /// \param contr Contravariant components of a vector.
    /// \param pos Current location.
    /// \return Covariant components of the transformed vector.
    template <class Coord>
    [[nodiscard]] auto contr_to_covar(Vector const &contr, Coord const &pos) const noexcept { return contr_to_covar(contr, self().covar_metric(pos)); }

    /// Vector transformation from covariant to contravariant components
    /// \tparam Coord Coordinate type.
    /// \param covar Covariant components of a vector.
    /// \param pos Current location.
    /// \return Contravariant components of the transformed vector.
    template <class Coord>
    [[nodiscard]] auto covar_to_contr(Vector const &covar, Coord const &pos) const noexcept { return covar_to_contr(covar, self().contr_metric(pos)); }

    /// Vector transformation from Cartesian to contravariant components
    /// \tparam Coord Coordinate type.
    /// \param cart Cartesian components of a vector.
    /// \param pos Current location.
    /// \return Contravariant components of the transformed vector.
    template <class Coord>
    [[nodiscard]] auto cart_to_contr(Vector const &cart, Coord const &pos) const noexcept { return cart_to_contr(cart, self().template contr_basis<0>(pos)); }

    /// Vector transformation from contravariant to Cartesian components
    /// \tparam Coord Coordinate type.
    /// \param contr Contravariant components of a vector.
    /// \param pos Current location.
    /// \return Cartesian components of the transformed vector.
    template <class Coord>
    [[nodiscard]] auto contr_to_cart(Vector const &contr, Coord const &pos) const noexcept { return contr_to_cart(contr, self().template covar_basis<0>(pos)); }

    /// Vector transformation from Cartesian to covaraint components
    /// \tparam Coord Coordinate type.
    /// \param cart Cartesian components of a vector.
    /// \param pos Current location.
    /// \return Covariant components of the transformed vector.
    template <class Coord>
    [[nodiscard]] auto cart_to_covar(Vector const &cart, Coord const &pos) const noexcept { return cart_to_covar(cart, self().template covar_basis<0>(pos)); }

    /// Vector transformation from covariant to Cartesian components
    /// \tparam Coord Coordinate type.
    /// \param covar Covariant components of a vector.
    /// \param pos Current location.
    /// \return Cartesian components of the transformed vector.
    template <class Coord>
    [[nodiscard]] auto covar_to_cart(Vector const &covar, Coord const &pos) const noexcept { return covar_to_cart(covar, self().template contr_basis<0>(pos)); }
};
} // namespace Detail
LIBPIC_NAMESPACE_END(1)
