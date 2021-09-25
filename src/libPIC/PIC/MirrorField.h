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
#include <PIC/Vector.h>

#include <cmath>

LIBPIC_BEGIN_NAMESPACE
namespace Detail {
template <class MirrorGeometry>
class MirrorField {
    [[nodiscard]] inline decltype(auto) self() const noexcept { return static_cast<MirrorGeometry const &>(*this); }

public:
    MirrorField() noexcept = default;

    [[nodiscard]] Vector B_div_B0(CartCoord const &pos) const noexcept
    {
        return { (1 + self().xi2() * pos.x * pos.x), 0, 0 };
    }
    [[nodiscard]] Vector B_div_B0(CartCoord const &pos, Real pos_y, Real pos_z) const noexcept
    {
        auto const xi2 = self().xi2();
        return {
            (1 + xi2 * pos.x * pos.x),
            (0 - xi2 * pos.x * pos_y),
            (0 - xi2 * pos.x * pos_z)
        };
    }

    [[nodiscard]] Vector B_div_B0(CurviCoord const &pos) const noexcept
    {
        auto const cos = std::cos(self().xi() * self().D1() * pos.q1);
        return { 1 / (cos * cos), 0, 0 };
    }
    [[nodiscard]] Vector B_div_B0(CurviCoord const &pos, Real pos_y, Real pos_z) const noexcept
    {
        return B_div_B0(self().cotrans(pos), pos_y, pos_z);
    }
};
} // namespace Detail
LIBPIC_END_NAMESPACE
