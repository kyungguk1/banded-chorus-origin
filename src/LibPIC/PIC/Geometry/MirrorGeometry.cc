/*
 * Copyright (c) 2021-2022, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "MirrorGeometry.h"
#include "../UTL/lippincott.h"

#include <cmath>
#include <limits>
#include <stdexcept>

LIBPIC_NAMESPACE_BEGIN(1)
namespace {
constexpr auto quiet_nan = std::numeric_limits<Real>::quiet_NaN();
}
Detail::MirrorGeometry::MirrorGeometry() noexcept
: m_D{ quiet_nan }
, m_inv_D{ quiet_nan }
, m_xi{ quiet_nan }
, m_xi2{ quiet_nan }
, m_sqrt_g{ quiet_nan }
, m_det_gij{ quiet_nan }
{
}
Detail::MirrorGeometry::MirrorGeometry(Real const xi, Vector const &D)
: m_homogeneous{ xi < inhomogeneity_xi_threshold }
{
    if (xi < 0)
        throw std::invalid_argument{ std::string{ __PRETTY_FUNCTION__ } + " - negative xi" };
    if (D.x <= 0)
        throw std::invalid_argument{ std::string{ __PRETTY_FUNCTION__ } + " - non-positive D1" };
    if (D.y <= 0)
        throw std::invalid_argument{ std::string{ __PRETTY_FUNCTION__ } + " - non-positive D2" };
    if (D.z <= 0)
        throw std::invalid_argument{ std::string{ __PRETTY_FUNCTION__ } + " - non-positive D3" };

    m_D       = D;
    m_inv_D   = 1 / m_D;
    m_xi      = xi;
    m_xi2     = xi * xi;
    m_sqrt_g  = m_D.x * m_D.y * m_D.z;
    m_det_gij = m_sqrt_g * m_sqrt_g;

    if (m_homogeneous) {
        m_cart_to_curvi = &MirrorGeometry::template cart_to_curvi<true>;
        m_curvi_to_cart = &MirrorGeometry::template curvi_to_cart<true>;
    } else {
        m_cart_to_curvi = &MirrorGeometry::template cart_to_curvi<false>;
        m_curvi_to_cart = &MirrorGeometry::template curvi_to_cart<false>;
    }
}

template <bool homogeneous>
auto Detail::MirrorGeometry::cart_to_curvi(CartCoord const &pos) const noexcept -> CurviCoord
{
    if constexpr (homogeneous) {
        return CurviCoord{ pos.x * inv_D1() };
    } else {
        return CurviCoord{ std::atan(xi() * pos.x) / (xi() * D1()) };
    }
}
template <bool homogeneous>
auto Detail::MirrorGeometry::curvi_to_cart(CurviCoord const &pos) const noexcept -> CartCoord
{
    if constexpr (homogeneous) {
        return CartCoord{ pos.q1 * D1() };
    } else {
        auto const xiD1q1 = xi() * D1() * pos.q1;
#if defined(DEBUG)
        if (std::abs(xiD1q1) >= M_PI_2)
            fatal_error("|xi*D1*q1| cannot be larger than pi/2");
#endif
        return CartCoord{ std::tan(xiD1q1) / xi() };
    }
}
LIBPIC_NAMESPACE_END(1)
