/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC/Config.h>
#include <PIC/Predefined.h>

#include <stdexcept>
#include <tuple>
#include <utility>

LIBPIC_BEGIN_NAMESPACE
/// Common parameters for all plasma populations
///
struct PlasmaDesc {
    Real Oc;                          //!< Cyclotron frequency.
    Real op;                          //!< Plasma frequency.
    long number_of_source_smoothings; //!< The number of source smoothings.

    /// Construct a plasma description
    /// \param Oc Cyclotron frequency of this plasma species/component.
    /// \param op Plasma frequency of this plasma species/component.
    /// \param n_smooths An optional argument to set the number of source smoothing. Default is 0.
    /// \throw Throws std::invalid_argument if either Oc == 0 or op <= 0.
    ///
    constexpr PlasmaDesc(Real Oc, Real op, unsigned n_smooths = {})
    : Oc{ Oc }, op{ op }, number_of_source_smoothings{ n_smooths }
    {
        if (this->Oc == 0)
            throw std::invalid_argument{ "Oc should not be zero" };
        if (this->op <= 0)
            throw std::invalid_argument{ "op should be positive" };
    }

protected:
    PlasmaDesc() noexcept = default;

    [[nodiscard]] friend constexpr auto serialize(PlasmaDesc const &desc) noexcept
    {
        return std::make_tuple(desc.Oc, desc.op);
    }
    [[nodiscard]] friend constexpr bool operator==(PlasmaDesc const &lhs, PlasmaDesc const &rhs) noexcept
    {
        return serialize(lhs) == serialize(rhs);
    }
};

/// Charge-neutralizing electron fluid descriptor.
///
struct eFluidDesc : public PlasmaDesc {
    Real beta;  //!< Electron beta.
    Real gamma; //!< Specific heat ratio, gamma.

    /// Construct charge-neutralizing, massless electron fluid description
    /// \param desc Common plasma description.
    /// \param beta Electron plasma beta. Default is 0.
    /// \param closure Polytropic index in equation of state. Default is adiabatic.
    explicit constexpr eFluidDesc(PlasmaDesc const &desc, Real beta = {}, Closure closure = adiabatic)
    : PlasmaDesc(desc), beta{ beta }, gamma(closure / 10)
    {
        gamma /= closure % 10;
        if (this->beta < 0)
            throw std::invalid_argument{ "beta should be non-negative" };
    }

private:
    [[nodiscard]] friend constexpr auto serialize(eFluidDesc const &desc) noexcept
    {
        PlasmaDesc const &base = desc;
        return std::tuple_cat(serialize(base), std::make_tuple(desc.beta, desc.gamma));
    }
};

/// Parameter set for a cold plasma population
///
struct ColdPlasmaDesc : public PlasmaDesc {
    Real Vd; //!< Equilibrium parallel drift speed.

    // the explicit qualifier is to prevent an accidental construction of an empty object
    //
    explicit ColdPlasmaDesc() noexcept = default;

    /// Construct a cold plasma description
    /// \param desc Common plasma description.
    /// \param Vd Equilibrium parallel drift speed. Default is 0.
    /// \throw Any exception thrown by PlasmaDesc.
    ///
    explicit constexpr ColdPlasmaDesc(PlasmaDesc const &desc, Real Vd = 0)
    : PlasmaDesc(desc), Vd{ Vd }
    {
    }

private:
    [[nodiscard]] friend constexpr auto serialize(ColdPlasmaDesc const &desc) noexcept
    {
        PlasmaDesc const &base = desc;
        return std::tuple_cat(serialize(base), std::make_tuple(desc.Vd));
    }
    [[nodiscard]] friend constexpr bool operator==(ColdPlasmaDesc const &lhs, ColdPlasmaDesc const &rhs) noexcept
    {
        return serialize(lhs) == serialize(rhs);
    }
};

/// Common parameters for all kinetic plasma populations
///
struct KineticPlasmaDesc : public PlasmaDesc {
    long           Nc;                //!< The number of simulation particles per cell.
    ShapeOrder     shape_order;       //!< The order of the shape function.
    ParticleScheme scheme;            //!< Full-f or delta-f scheme.
    Real           initial_weight;    //!< Initial particle's delta-f weight.
    Real           marker_temp_ratio; //!< Relative fraction of marker particle temperature.

    // the explicit qualifier is to prevent an accidental construction of an empty object
    //
    explicit KineticPlasmaDesc() noexcept = default;

    /// Construct a kinetic plasma description
    /// \param desc Common plasma description.
    /// \param Nc Number of simulation particles.
    /// \param shape_order Simulation particle shape order.
    /// \param scheme Whether to evolve full or delta VDF. Default is full_f.
    /// \param initial_weight Initial weight of delta-f particles. Default is 0.
    /// \param marker_temp_ratio Relative fraction of marker particle's temperature.
    ///                          Must be positive. Default is 1.
    /// \throw Any exception thrown by PlasmaDesc, and if Nc == 0.
    ///
    constexpr KineticPlasmaDesc(PlasmaDesc const &desc, unsigned Nc, ShapeOrder shape_order,
                                ParticleScheme scheme = full_f, Real initial_weight = 0, Real marker_temp_ratio = 1)
    : PlasmaDesc(desc)
    , Nc{ Nc }
    , shape_order{ shape_order }
    , scheme{ scheme }
    , initial_weight{ full_f == scheme ? 0 : initial_weight }
    , marker_temp_ratio{ full_f == scheme ? 1 : marker_temp_ratio }
    {
        if (this->Nc <= 0)
            throw std::invalid_argument{ "Nc should be positive" };
        if (this->initial_weight < 0 || this->initial_weight >= 1)
            throw std::invalid_argument{ "initial weight should be between 0 and 1 (exclusive)" };
        if (this->marker_temp_ratio <= 0)
            throw std::invalid_argument{ "relative fraction of marker particle's temperature must be a positive number" };
    }

private:
    [[nodiscard]] friend constexpr auto serialize(KineticPlasmaDesc const &desc) noexcept
    {
        PlasmaDesc const &base = desc;
        return std::tuple_cat(serialize(base),
                              std::make_tuple(desc.Nc, desc.scheme, desc.initial_weight, desc.marker_temp_ratio));
    }
    [[nodiscard]] friend constexpr bool operator==(KineticPlasmaDesc const &lhs, KineticPlasmaDesc const &rhs) noexcept
    {
        return serialize(lhs) == serialize(rhs);
    }
};

/// Parameters for a bi-Maxwellian plasma population
///
struct BiMaxPlasmaDesc : public KineticPlasmaDesc {
    Real beta1; //!< The parallel component of plasma beta.
    Real T2_T1; //!< The ratio of the perpendicular to parallel temperatures.
    Real Vd;    //!< Equilibrium parallel drift speed.

    /// Construct a bi-Maxwellian plasma description.
    /// \param desc Kinetic plasma description.
    /// \param beta1 Parallel plasma beta.
    /// \param T2_T1 Perpendicular-to-parallel temperature ratio. Default is 1.
    /// \param Vd Parallel bulk drift velocity. Default is 0.
    /// \throw Any exception thrown by KineticPlasmaDesc, and if either beta1 <= 0 or T2_T1 <= 0.
    ///
    constexpr BiMaxPlasmaDesc(KineticPlasmaDesc const &desc, Real beta1, Real T2_T1 = 1, Real Vd = 0)
    : KineticPlasmaDesc(desc), beta1{ beta1 }, T2_T1{ T2_T1 }, Vd{ Vd }
    {
        if (this->beta1 <= 0)
            throw std::invalid_argument{ "beta1 should be positive" };
        if (this->T2_T1 <= 0)
            throw std::invalid_argument{ "T2_T1 should be positive" };
    }

private:
    [[nodiscard]] friend constexpr auto serialize(BiMaxPlasmaDesc const &desc) noexcept
    {
        KineticPlasmaDesc const &base = desc;
        return std::tuple_cat(serialize(base), std::make_tuple(desc.beta1, desc.T2_T1, desc.Vd));
    }
    [[nodiscard]] friend constexpr bool operator==(BiMaxPlasmaDesc const &lhs, BiMaxPlasmaDesc const &rhs) noexcept
    {
        return serialize(lhs) == serialize(rhs);
    }
};

/// Parameters for a loss-cone distribution plasma population
/// \details The perpendicular component of the loss-cone is given by
///          f_perp = ((1 - Δβ)*exp(-x^2) - (1 - Δ)*exp(-x^2/β)) / (1 - β)*π*θ2^2
/// where x = v2/θ2.
/// The effective perpendicular temperature is 2*T2 = 1 + (1 - Δ)*β.
///
struct LossconePlasmaDesc : public BiMaxPlasmaDesc {
    Real Delta; // Loss-cone VDF Δ parameter.
    Real beta;  // Loss-cone VDF β parameter.

    /// Construct a loss-cone plasma description
    /// \details In this version, the effective temperatures are used to derive the necessary
    /// parameters.
    /// \param desc A bi-Maxwellian plasma description.
    /// \param Delta ∆ parameter. Default is 1, in which case it becomes a bi-Maxwellian.
    /// \param beta β parameter. Default is 1.
    /// \throw Any exception thrown by BiMaxPlasmaDesc, and
    ///        if either ∆ lies outside the range [0, 1] or β <= 0.
    ///
    explicit constexpr LossconePlasmaDesc(BiMaxPlasmaDesc const &desc, Real Delta = 1, Real beta = 1)
    : BiMaxPlasmaDesc(desc), Delta{ Delta }, beta{ beta }
    {
        if (this->Delta < 0 || this->Delta > 1)
            throw std::invalid_argument{ "Losscone.Delta should be in the range of [0, 1]" };
        if (this->beta <= 0)
            throw std::invalid_argument{ "Losscone.beta should be positive" };
    }

    /// Construct a loss-cone plasma description
    /// \details In this version, the necessary parameters are explicitly specified.
    /// \param desc A kinetic plasma description.
    /// \param beta1 Parallel plasma beta.
    /// \param vth_ratio A positive number for the ratio θ2^2/θ1^2. Default is 1.
    /// \param Db A pair of {∆, β}. Default is {1, 1}.
    /// \param Vd Parallel bulk drift velocity. Default is 0.
    /// \throw Any exception thrown by BiMaxPlasmaDesc, and
    ///        if either ∆ lies outside the range [0, 1] or β <= 0.
    ///
    explicit constexpr LossconePlasmaDesc(KineticPlasmaDesc const &desc, Real beta1, Real vth_ratio = 1,
                                          std::pair<Real, Real> Db = { 1, 1 }, Real Vd = 0)
    : LossconePlasmaDesc({ desc, beta1, (1 + (1 - Db.first) * Db.second) * vth_ratio, Vd }, Db.first, Db.second)
    {
    }

private:
    [[nodiscard]] friend constexpr auto serialize(LossconePlasmaDesc const &desc) noexcept
    {
        BiMaxPlasmaDesc const &base = desc;
        return std::tuple_cat(serialize(base), std::make_tuple(desc.Delta, desc.beta));
    }
    [[nodiscard]] friend constexpr bool operator==(LossconePlasmaDesc const &lhs, LossconePlasmaDesc const &rhs) noexcept
    {
        return serialize(lhs) == serialize(rhs);
    }
};
LIBPIC_END_NAMESPACE
