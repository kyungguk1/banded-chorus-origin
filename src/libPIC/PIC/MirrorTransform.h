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
#include <PIC/Tensor.h>
#include <PIC/Vector.h>

LIBPIC_BEGIN_NAMESPACE
namespace Detail {
template <class MirrorGeometry>
class MirrorTransform {
    [[nodiscard]] inline decltype(auto) self() const noexcept { return static_cast<MirrorGeometry const &>(*this); }

protected:
    MirrorTransform() noexcept = default;

public:
    [[nodiscard]] static Vector contr_to_covar(Vector const &contr, Tensor const &covar_metric) noexcept { return dot(contr, covar_metric); }
    [[nodiscard]] static Vector covar_to_contr(Vector const &covar, Tensor const &contr_metric) noexcept { return dot(covar, contr_metric); }
    [[nodiscard]] static Vector cart_to_contr(Vector const &cart, Tensor const &contr_bases) noexcept { return dot(contr_bases, cart); }
    [[nodiscard]] static Vector contr_to_cart(Vector const &contr, Tensor const &covar_bases) noexcept { return dot(contr, covar_bases); }
    [[nodiscard]] static Vector cart_to_covar(Vector const &cart, Tensor const &covar_bases) noexcept { return dot(covar_bases, cart); }
    [[nodiscard]] static Vector covar_to_cart(Vector const &covar, Tensor const &contr_bases) noexcept { return dot(covar, contr_bases); }

    template <class Coord>
    [[nodiscard]] auto contr_to_covar(Vector const &contr, Coord const &pos) const noexcept { return contr_to_covar(contr, self().covar_metric(pos)); }
    template <class Coord>
    [[nodiscard]] auto covar_to_contr(Vector const &covar, Coord const &pos) const noexcept { return covar_to_contr(covar, self().contr_metric(pos)); }

    template <class Coord>
    [[nodiscard]] auto cart_to_contr(Vector const &cart, Coord const &pos) const noexcept { return cart_to_contr(cart, self().template contr_basis<0>(pos)); }
    template <class Coord>
    [[nodiscard]] auto contr_to_cart(Vector const &contr, Coord const &pos) const noexcept { return contr_to_cart(contr, self().template covar_basis<0>(pos)); }

    template <class Coord>
    [[nodiscard]] auto cart_to_covar(Vector const &cart, Coord const &pos) const noexcept { return cart_to_covar(cart, self().template covar_basis<0>(pos)); }
    template <class Coord>
    [[nodiscard]] auto covar_to_cart(Vector const &covar, Coord const &pos) const noexcept { return covar_to_cart(covar, self().template contr_basis<0>(pos)); }
};
} // namespace Detail
LIBPIC_END_NAMESPACE
