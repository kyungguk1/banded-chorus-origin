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

LIBPIC_BEGIN_NAMESPACE
class MirrorGeometry;

namespace Detail {
class MirrorCotrans {
    [[nodiscard]] inline decltype(auto) self() const noexcept;

protected:
    MirrorCotrans() noexcept = default;
    explicit MirrorCotrans(bool homogeneous) noexcept;

public:
    [[nodiscard]] CurviCoord cotrans(CartCoord const &cart) const noexcept { return (this->*m_cart_to_curvi)(cart); };
    [[nodiscard]] CartCoord  cotrans(CurviCoord const &curvi) const noexcept { return (this->*m_curvi_to_cart)(curvi); };

private:
    template <bool homogeneous>
    auto cart_to_curvi(CartCoord const &) const noexcept;
    template <bool homogeneous>
    auto curvi_to_cart(CurviCoord const &) const noexcept;

    CurviCoord (MirrorCotrans::*m_cart_to_curvi)(CartCoord const &) const noexcept;
    CartCoord (MirrorCotrans::*m_curvi_to_cart)(CurviCoord const &) const noexcept;
};
} // namespace Detail
LIBPIC_END_NAMESPACE
