/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC/RelativisticLossconeVDF.h>
#include <PIC/RelativisticMaxwellianVDF.h>
#include <PIC/RelativisticPartialShellVDF.h>
#include <PIC/RelativisticTestParticleVDF.h>

#include <stdexcept>
#include <type_traits>
#include <utility>
#include <variant>

LIBPIC_BEGIN_NAMESPACE
class RelativisticVDFVariant {
    // visitor overload facility
    //
    template <class Ret, class Lambda>
    struct Overload : Lambda {
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
    using variant_t = std::variant<std::monostate, RelativisticMaxwellianVDF, RelativisticLossconeVDF, RelativisticTestParticleVDF, RelativisticPartialShellVDF>;
    using Particle  = RelativisticParticle;

    // ctor's
    //
    RelativisticVDFVariant() = default;

    template <class... Args>
    [[nodiscard]] static RelativisticVDFVariant make(BiMaxPlasmaDesc const &desc, Args &&...args) noexcept(
        std::is_nothrow_constructible_v<RelativisticMaxwellianVDF, decltype(desc), Args...>)
    {
        static_assert(std::is_constructible_v<RelativisticMaxwellianVDF, decltype(desc), Args...>);
        return { std::in_place_type<RelativisticMaxwellianVDF>, desc, std::forward<Args>(args)... };
    }
    template <class... Args>
    [[nodiscard]] static RelativisticVDFVariant make(LossconePlasmaDesc const &desc, Args &&...args) noexcept(
        std::is_nothrow_constructible_v<RelativisticLossconeVDF, decltype(desc), Args...>)
    {
        static_assert(std::is_constructible_v<RelativisticLossconeVDF, decltype(desc), Args...>);
        return { std::in_place_type<RelativisticLossconeVDF>, desc, std::forward<Args>(args)... };
    }
    template <unsigned N, class... Args>
    [[nodiscard]] static RelativisticVDFVariant make(TestParticleDesc<N> const &desc, Args &&...args) noexcept(
        std::is_nothrow_constructible_v<RelativisticTestParticleVDF, decltype(desc), Args...>)
    {
        static_assert(std::is_constructible_v<RelativisticTestParticleVDF, decltype(desc), Args...>);
        return { std::in_place_type<RelativisticTestParticleVDF>, desc, std::forward<Args>(args)... };
    }
    template <class... Args>
    [[nodiscard]] static RelativisticVDFVariant make(PartialShellPlasmaDesc const &desc, Args &&...args) noexcept(
        std::is_nothrow_constructible_v<RelativisticPartialShellVDF, decltype(desc), Args...>)
    {
        static_assert(std::is_constructible_v<RelativisticPartialShellVDF, decltype(desc), Args...>);
        return { std::in_place_type<RelativisticPartialShellVDF>, desc, std::forward<Args>(args)... };
    }

    template <class... Args>
    [[nodiscard]] decltype(auto) emplace(BiMaxPlasmaDesc const &desc, Args &&...args) noexcept(
        std::is_nothrow_constructible_v<RelativisticMaxwellianVDF, decltype(desc), Args...>)
    {
        static_assert(std::is_constructible_v<RelativisticMaxwellianVDF, decltype(desc), Args...>);
        return var.emplace<RelativisticMaxwellianVDF>(desc, std::forward<Args>(args)...);
    }
    template <class... Args>
    [[nodiscard]] decltype(auto) emplace(LossconePlasmaDesc const &desc, Args &&...args) noexcept(
        std::is_nothrow_constructible_v<RelativisticLossconeVDF, decltype(desc), Args...>)
    {
        static_assert(std::is_constructible_v<RelativisticLossconeVDF, decltype(desc), Args...>);
        return var.emplace<RelativisticLossconeVDF>(desc, std::forward<Args>(args)...);
    }
    template <unsigned N, class... Args>
    [[nodiscard]] decltype(auto) emplace(TestParticleDesc<N> const &desc, Args &&...args) noexcept(
        std::is_nothrow_constructible_v<RelativisticTestParticleVDF, decltype(desc), Args...>)
    {
        static_assert(std::is_constructible_v<RelativisticTestParticleVDF, decltype(desc), Args...>);
        return var.emplace<RelativisticTestParticleVDF>(desc, std::forward<Args>(args)...);
    }
    template <class... Args>
    [[nodiscard]] decltype(auto) emplace(PartialShellPlasmaDesc const &desc, Args &&...args) noexcept(
        std::is_nothrow_constructible_v<RelativisticPartialShellVDF, decltype(desc), Args...>)
    {
        static_assert(std::is_constructible_v<RelativisticPartialShellVDF, decltype(desc), Args...>);
        return var.emplace<RelativisticPartialShellVDF>(desc, std::forward<Args>(args)...);
    }

    // method dispatch
    //
    [[nodiscard]] KineticPlasmaDesc const &plasma_desc() const noexcept
    {
        using Ret          = decltype(plasma_desc());
        constexpr auto vis = make_vis<Ret>([](auto const &alt) -> Ret {
            return alt.plasma_desc();
        });
        return std::visit(vis, var);
    }
    [[nodiscard]] Particle emit() const
    {
        using Ret          = decltype(emit());
        constexpr auto vis = make_vis<Ret>([](auto const &alt) -> Ret {
            return alt.emit();
        });
        return std::visit(vis, var);
    }
    [[nodiscard]] std::vector<Particle> emit(unsigned long n) const
    {
        using Ret      = decltype(emit(n));
        const auto vis = make_vis<Ret>([n](auto const &alt) -> Ret {
            return alt.emit(n);
        });
        return std::visit(vis, var);
    }

    [[nodiscard]] Scalar n0(CurviCoord const &pos) const
    {
        using Ret      = decltype(n0(pos));
        const auto vis = make_vis<Ret>([&pos](auto const &alt) -> Ret {
            return alt.n0(pos);
        });
        return std::visit(vis, var);
    }
    [[nodiscard]] Vector nV0(CurviCoord const &pos) const
    {
        using Ret      = decltype(nV0(pos));
        const auto vis = make_vis<Ret>([&pos](auto const &alt) -> Ret {
            return alt.nV0(pos);
        });
        return std::visit(vis, var);
    }
    [[nodiscard]] FourTensor nuv0(CurviCoord const &pos) const
    {
        using Ret      = decltype(nuv0(pos));
        const auto vis = make_vis<Ret>([&pos](auto const &alt) -> Ret {
            return alt.nuv0(pos);
        });
        return std::visit(vis, var);
    }

    [[nodiscard]] Real weight(Particle const &ptl) const
    {
        using Ret      = decltype(weight(ptl));
        const auto vis = make_vis<Ret>([&ptl](auto const &alt) -> Ret {
            return alt.weight(ptl);
        });
        return std::visit(vis, var);
    }

    [[nodiscard]] Real Nrefcell_div_Ntotal() const
    {
        using Ret      = decltype(Nrefcell_div_Ntotal());
        const auto vis = make_vis<Ret>([](auto const &alt) -> Ret {
            return alt.Nrefcell_div_Ntotal();
        });
        return std::visit(vis, var);
    }

private:
    template <class VDF, class... Args>
    RelativisticVDFVariant(std::in_place_type_t<VDF> type, Args &&...args) noexcept(std::is_nothrow_constructible_v<VDF, Args...>)
    : var{ type, std::forward<Args>(args)... }
    {
    }

    variant_t var;
};
LIBPIC_END_NAMESPACE
