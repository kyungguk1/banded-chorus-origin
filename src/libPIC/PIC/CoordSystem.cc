/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "CoordSystem.h"
#include <PIC/lippincott.h>
#include <cmath>
#include <stdexcept>

LIBPIC_BEGIN_NAMESPACE
CoordSystem::CoordSystem(Real const xi, Real const D1)
{
    if (xi < 0)
        throw std::invalid_argument{ std::string{ __PRETTY_FUNCTION__ } + " - negative xi" };
    if (D1 <= 0)
        throw std::invalid_argument{ std::string{ __PRETTY_FUNCTION__ } + " - non-positive D1" };

    m_D           = { D1, 1, 1 };
    m_xi          = xi;
    m_xi2         = xi * xi;
    m_sqrt_g      = m_D.x * m_D.y * m_D.z;
    m_det_gij     = m_sqrt_g * m_sqrt_g;
    m_homogeneous = m_xi2 < 1e-10;

    // dispatch
    if (m_homogeneous) {
        m_cart_to_curvi = &CoordSystem::template cart_to_curvi<true>;
        m_curvi_to_cart = &CoordSystem::template curvi_to_cart<true>;
    } else {
        m_cart_to_curvi = &CoordSystem::template cart_to_curvi<false>;
        m_curvi_to_cart = &CoordSystem::template curvi_to_cart<false>;
    }
}

template <bool homogeneous>
auto CoordSystem::cart_to_curvi(CartCoord const &pos) const noexcept -> CurviCoord
{
    if constexpr (homogeneous) {
        return CurviCoord{ pos.x / D1() };
    } else {
        return CurviCoord{ std::atan(xi() * pos.x) / (xi() * D1()) };
    }
}
template <bool homogeneous>
auto CoordSystem::curvi_to_cart(CurviCoord const &pos) const noexcept -> CartCoord
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

LIBPIC_END_NAMESPACE
