/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC/LossconeVDF.h>
#include <PIC/MaxwellianVDF.h>

#include <stdexcept>
#include <type_traits>
#include <utility>
#include <variant>

LIBPIC_BEGIN_NAMESPACE
class VDFVariant {
    // visitor overload facility
    //
    template <class Ret, class Lambda> struct Overload : Lambda {
        using Lambda::operator();

        Ret operator()(std::monostate const &) const
        {
            throw std::domain_error{ __PRETTY_FUNCTION__ };
        }
    };

    template <class Ret, class Vis>
    [[nodiscard]] static constexpr Overload<Ret, Vis> make_vis(Vis vis)
    {
        return { std::move(vis) };
    }

public:
    using variant_t = std::variant<std::monostate, MaxwellianVDF, LossconeVDF>;

    // ctor's
    //
    VDFVariant() = default;

    template <class VDF, class... Args>
    [[nodiscard]] static VDFVariant
    create(Args &&...args) noexcept(std::is_nothrow_constructible_v<VDF, Args...>)
    {
        static_assert(std::is_constructible_v<VDF, Args...>, "VDF object is not constructible");
        return { std::in_place_type<VDF>, std::forward<Args>(args)... };
    }

    template <class VDF, class... Args>
    [[nodiscard]] decltype(auto)
    emplace(Args &&...args) noexcept(std::is_nothrow_constructible_v<VDF, Args...>)
    {
        static_assert(std::is_constructible_v<VDF, Args...>, "VDF object is not constructible");
        return var.emplace<VDF>(std::forward<Args>(args)...);
    }

    // method dispatch
    //
    [[nodiscard]] Particle emit() const
    {
        constexpr auto vis = make_vis<decltype(emit())>([](auto const &alt) {
            return alt.emit();
        });
        return std::visit(vis, var);
    }
    [[nodiscard]] std::vector<Particle> emit(unsigned n) const
    {
        const auto vis = make_vis<decltype(emit(n))>([n](auto const &alt) {
            return alt.emit(n);
        });
        return std::visit(vis, var);
    }

    [[nodiscard]] Scalar n0(Real pos_x) const
    {
        const auto vis = make_vis<decltype(n0(pos_x))>([pos_x](auto const &alt) {
            return alt.n0(pos_x);
        });
        return std::visit(vis, var);
    }
    [[nodiscard]] Vector nV0(Real pos_x) const
    {
        const auto vis = make_vis<decltype(nV0(pos_x))>([pos_x](auto const &alt) {
            return alt.nV0(pos_x);
        });
        return std::visit(vis, var);
    }
    [[nodiscard]] Tensor nvv0(Real pos_x) const
    {
        const auto vis = make_vis<decltype(nvv0(pos_x))>([pos_x](auto const &alt) {
            return alt.nvv0(pos_x);
        });
        return std::visit(vis, var);
    }

    [[nodiscard]] Real delta_f(Particle const &ptl) const
    {
        const auto vis = make_vis<decltype(delta_f(ptl))>([&ptl](auto const &alt) {
            return alt.delta_f(ptl);
        });
        return std::visit(vis, var);
    }
    [[nodiscard]] Real weight(Particle const &ptl) const
    {
        const auto vis = make_vis<decltype(weight(ptl))>([&ptl](auto const &alt) {
            return alt.weight(ptl);
        });
        return std::visit(vis, var);
    }

private:
    template <class VDF, class... Args>
    VDFVariant(std::in_place_type_t<VDF> type,
               Args &&...args) noexcept(std::is_nothrow_constructible_v<VDF, Args...>)
    : var{ type, std::forward<Args>(args)... }
    {
    }

    variant_t var;
};
LIBPIC_END_NAMESPACE
