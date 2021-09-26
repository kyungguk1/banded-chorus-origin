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

#include <cmath>
#include <type_traits>

LIBPIC_BEGIN_NAMESPACE
namespace Detail {
template <class MirrorGeometry>
class MirrorBasis {
    [[nodiscard]] inline decltype(auto) self() const noexcept { return static_cast<MirrorGeometry const &>(*this); }

    [[nodiscard]] inline static auto pow2(Real const x) noexcept { return x * x; }

    [[nodiscard]] Vector covar_basis(CartCoord const &pos, std::integral_constant<long, 1>) const noexcept { return { self().D1() * (1 + self().xi2() * pow2(pos.x)), 0, 0 }; }
    [[nodiscard]] Vector covar_basis(CartCoord const &pos, std::integral_constant<long, 2>) const noexcept { return { 0, self().D2() / std::sqrt(1 + self().xi2() * pow2(pos.x)), 0 }; }
    [[nodiscard]] Vector covar_basis(CartCoord const &pos, std::integral_constant<long, 3>) const noexcept { return { 0, 0, self().D3() / std::sqrt(1 + self().xi2() * pow2(pos.x)) }; }
    [[nodiscard]] Tensor covar_basis(CartCoord const &pos, std::integral_constant<long, 0>) const noexcept
    {
        auto const tmp      = 1 + self().xi2() * pow2(pos.x);
        auto const sqrt_tmp = std::sqrt(tmp);
        return { self().D1() * tmp, self().D2() / sqrt_tmp, self().D3() / sqrt_tmp, 0, 0, 0 };
    }

    [[nodiscard]] Vector covar_basis(CurviCoord const &pos, std::integral_constant<long, 1>) const noexcept { return { self().D1() / pow2(std::cos(self().xi() * self().D1() * pos.q1)), 0, 0 }; }
    [[nodiscard]] Vector covar_basis(CurviCoord const &pos, std::integral_constant<long, 2>) const noexcept { return { 0, self().D2() * std::cos(self().xi() * self().D1() * pos.q1), 0 }; }
    [[nodiscard]] Vector covar_basis(CurviCoord const &pos, std::integral_constant<long, 3>) const noexcept { return { 0, 0, self().D3() * std::cos(self().xi() * self().D1() * pos.q1) }; }
    [[nodiscard]] Tensor covar_basis(CurviCoord const &pos, std::integral_constant<long, 0>) const noexcept
    {
        auto const tmp = std::cos(self().xi() * self().D1() * pos.q1);
        return { self().D1() / pow2(tmp), self().D2() * tmp, self().D3() * tmp, 0, 0, 0 };
    }

    [[nodiscard]] Vector contr_basis(CartCoord const &pos, std::integral_constant<long, 1>) const noexcept { return { self().inv_D1() / (1 + self().xi2() * pow2(pos.x)), 0, 0 }; }
    [[nodiscard]] Vector contr_basis(CartCoord const &pos, std::integral_constant<long, 2>) const noexcept { return { 0, std::sqrt(1 + self().xi2() * pow2(pos.x)) * self().inv_D2(), 0 }; }
    [[nodiscard]] Vector contr_basis(CartCoord const &pos, std::integral_constant<long, 3>) const noexcept { return { 0, 0, std::sqrt(1 + self().xi2() * pow2(pos.x)) * self().inv_D3() }; }
    [[nodiscard]] Tensor contr_basis(CartCoord const &pos, std::integral_constant<long, 0>) const noexcept
    {
        auto const tmp      = 1 + self().xi2() * pow2(pos.x);
        auto const sqrt_tmp = std::sqrt(tmp);
        return { self().inv_D1() / tmp, sqrt_tmp * self().inv_D2(), sqrt_tmp * self().inv_D3(), 0, 0, 0 };
    }

    [[nodiscard]] Vector contr_basis(CurviCoord const &pos, std::integral_constant<long, 1>) const noexcept { return { self().inv_D1() * pow2(std::cos(self().xi() * self().D1() * pos.q1)), 0, 0 }; }
    [[nodiscard]] Vector contr_basis(CurviCoord const &pos, std::integral_constant<long, 2>) const noexcept { return { 0, self().inv_D2() / std::cos(self().xi() * self().D1() * pos.q1), 0 }; }
    [[nodiscard]] Vector contr_basis(CurviCoord const &pos, std::integral_constant<long, 3>) const noexcept { return { 0, 0, self().inv_D3() / std::cos(self().xi() * self().D1() * pos.q1) }; }
    [[nodiscard]] Tensor contr_basis(CurviCoord const &pos, std::integral_constant<long, 0>) const noexcept
    {
        auto const tmp = std::cos(self().xi() * self().D1() * pos.q1);
        return { self().inv_D1() * pow2(tmp), self().inv_D2() / tmp, self().inv_D3() / tmp, 0, 0, 0 };
    }

    template <long i, class Coord>
    [[nodiscard]] static constexpr auto mfa_basis(Coord const &, std::integral_constant<long, i>) noexcept
    {
        if constexpr (i == 1)
            return Vector{ 1, 0, 0 };
        if constexpr (i == 2)
            return Vector{ 0, 1, 0 };
        if constexpr (i == 3)
            return Vector{ 0, 0, 1 };
        if constexpr (i == 0)
            return Tensor{ 1, 1, 1, 0, 0, 0 };
    }

protected:
    MirrorBasis() noexcept = default;

public:
    [[nodiscard]] Tensor covar_metric(CartCoord const &pos) const noexcept
    {
        auto const tmp = 1 + self().xi2() * pow2(pos.x);
        return { pow2(self().D1() * tmp), pow2(self().D2()) / tmp, pow2(self().D3()) / tmp, 0, 0, 0 };
    }
    [[nodiscard]] Tensor covar_metric(CurviCoord const &pos) const noexcept
    {
        auto const tmp = pow2(std::cos(self().xi() * self().D1() * pos.q1));
        return { pow2(self().D1() / tmp), pow2(self().D2()) * tmp, pow2(self().D3()) * tmp, 0, 0, 0 };
    }

    [[nodiscard]] Tensor contr_metric(CartCoord const &pos) const noexcept
    {
        auto const tmp = 1 + self().xi2() * pow2(pos.x);
        return { pow2(self().inv_D1() / tmp), tmp * pow2(self().inv_D2()), tmp * pow2(self().inv_D3()), 0, 0, 0 };
    }
    [[nodiscard]] Tensor contr_metric(CurviCoord const &pos) const noexcept
    {
        auto const tmp = pow2(std::cos(self().xi() * self().D1() * pos.q1));
        return { pow2(self().inv_D1() * tmp), pow2(self().inv_D2()) / tmp, pow2(self().inv_D3()) / tmp, 0, 0, 0 };
    }

    template <long i>
    [[nodiscard]] auto covar_basis(CartCoord const &pos) const noexcept
    {
        static_assert(i >= 0 && i <= 3, "invalid index range");
        return covar_basis(pos, std::integral_constant<long, i>{});
    }
    template <long i>
    [[nodiscard]] auto covar_basis(CurviCoord const &pos) const noexcept
    {
        static_assert(i >= 0 && i <= 3, "invalid index range");
        return covar_basis(pos, std::integral_constant<long, i>{});
    }

    template <long i>
    [[nodiscard]] auto contr_basis(CartCoord const &pos) const noexcept
    {
        static_assert(i >= 0 && i <= 3, "invalid index range");
        return contr_basis(pos, std::integral_constant<long, i>{});
    }
    template <long i>
    [[nodiscard]] auto contr_basis(CurviCoord const &pos) const noexcept
    {
        static_assert(i >= 0 && i <= 3, "invalid index range");
        return contr_basis(pos, std::integral_constant<long, i>{});
    }

    template <long i>
    [[nodiscard]] static constexpr auto mfa_basis(CartCoord const &pos) noexcept
    {
        static_assert(i >= 0 && i <= 3, "invalid index range");
        return mfa_basis(pos, std::integral_constant<long, i>{});
    }
    template <long i>
    [[nodiscard]] static constexpr auto mfa_basis(CurviCoord const &pos) noexcept
    {
        static_assert(i >= 0 && i <= 3, "invalid index range");
        return mfa_basis(pos, std::integral_constant<long, i>{});
    }
};
} // namespace Detail
LIBPIC_END_NAMESPACE
