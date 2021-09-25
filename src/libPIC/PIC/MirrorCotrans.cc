/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "MirrorCotrans.h"
#include "MirrorGeometry.h"
#include "lippincott.h"
#include <cmath>

LIBPIC_BEGIN_NAMESPACE
auto Detail::MirrorCotrans::self() const noexcept -> MirrorGeometry const &
{
    return static_cast<MirrorGeometry const &>(*this);
}

template <bool homogeneous>
auto Detail::MirrorCotrans::cart_to_curvi(CartCoord const &pos) const noexcept
{
    if constexpr (homogeneous) {
        return CurviCoord{ pos.x / self().D1() };
    } else {
        return CurviCoord{ std::atan(self().xi() * pos.x) / (self().xi() * self().D1()) };
    }
}
template <bool homogeneous>
auto Detail::MirrorCotrans::curvi_to_cart(CurviCoord const &pos) const noexcept
{
    if constexpr (homogeneous) {
        return CartCoord{ pos.q1 * self().D1() };
    } else {
        auto const xiD1q1 = self().xi() * self().D1() * pos.q1;
#if defined(DEBUG)
        if (std::abs(xiD1q1) >= M_PI_2)
            fatal_error("|xi*D1*q1| cannot be larger than pi/2");
#endif
        return CartCoord{ std::tan(xiD1q1) / self().xi() };
    }
}
Detail::MirrorCotrans::MirrorCotrans(bool homogeneous) noexcept
{
    if (homogeneous) {
        m_cart_to_curvi = &MirrorCotrans::template cart_to_curvi<true>;
        m_curvi_to_cart = &MirrorCotrans::template curvi_to_cart<true>;
    } else {
        m_cart_to_curvi = &MirrorCotrans::template cart_to_curvi<false>;
        m_curvi_to_cart = &MirrorCotrans::template curvi_to_cart<false>;
    }
}
LIBPIC_END_NAMESPACE
