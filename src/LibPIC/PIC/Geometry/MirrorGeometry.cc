/*
 * Copyright (c) 2021-2022, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "MirrorGeometry.h"
#include "../UTL/lippincott.h"

#include <cmath>

LIBPIC_NAMESPACE_BEGIN(1)
Detail::MirrorGeometry::MirrorGeometry(Real const xi, Vector const &D)
: CurviBasis{ xi, D }
, MFABasis{ xi, D }
, m_D{ D }
, m_xi{ xi }
, m_sqrt_g{ m_D.x * m_D.y * m_D.z }
, m_det_gij{ m_sqrt_g * m_sqrt_g }
, m_homogeneous{ xi < inhomogeneity_xi_threshold }
{
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
    if constexpr (homogeneous)
        return CurviCoord{ pos.x / D1() };
    else
        return CurviCoord{ std::atan(xi() * pos.x) / (xi() * D1()) };
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
