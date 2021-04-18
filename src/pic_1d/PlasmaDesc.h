//
// Copyright (c) 2019, Kyungguk Min
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// - Redistributions of source code must retain the above copyright notice, this
//  list of conditions and the following disclaimer.
//
// - Redistributions in binary form must reproduce the above copyright notice,
//  this list of conditions and the following disclaimer in the documentation
//  and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#ifndef PlasmaDesc_h
#define PlasmaDesc_h

#include "./Macros.h"
#include "./Predefined.h"

#include <stdexcept>
#include <tuple>
#include <utility>

PIC1D_BEGIN_NAMESPACE
/// Common parameters for all plasmas.
///
struct [[nodiscard]] PlasmaDesc {
    Real Oc;                          //!< Cyclotron frequency.
    Real op;                          //!< Plasma frequency.
    long number_of_source_smoothings; //!< The number of source smoothings.
    //
    constexpr PlasmaDesc(Real Oc, Real op, unsigned n_smooths = {})
    : Oc{Oc}, op{op}, number_of_source_smoothings{n_smooths}
    {
        if (this->Oc == 0)
            throw std::invalid_argument{"Oc should not be zero"};
        if (this->op <= 0)
            throw std::invalid_argument{"op should be positive"};
    }

protected:
    explicit PlasmaDesc() noexcept = default;

private:
    [[nodiscard]] friend constexpr auto serialize(PlasmaDesc const &desc) noexcept
    {
        return std::make_tuple(desc.Oc, desc.op);
    }
};

/// Cold plasma descriptor.
///
struct [[nodiscard]] ColdPlasmaDesc : public PlasmaDesc {
    Real Vd; //!< Equilibrium parallel drift speed.
    //
    explicit ColdPlasmaDesc() noexcept = default;
    constexpr ColdPlasmaDesc(PlasmaDesc const &desc, Real Vd) : PlasmaDesc(desc), Vd{Vd} {}
    explicit constexpr ColdPlasmaDesc(PlasmaDesc const &desc) : ColdPlasmaDesc(desc, {}) {}

private:
    [[nodiscard]] friend constexpr auto serialize(ColdPlasmaDesc const &desc) noexcept
    {
        PlasmaDesc const &base = desc;
        return std::tuple_cat(serialize(base), std::make_tuple(desc.Vd));
    }
};

//
// MARK:- Kinetic Plasmas
//

/// Common parameters for all kinetic plasmas.
///
struct [[nodiscard]] KineticPlasmaDesc : public PlasmaDesc {
    long           Nc;          //!< The number of simulation particles per cell.
    ShapeOrder     shape_order; //!< The order of the shape function.
    ParticleScheme scheme;      //!< Full-f or delta-f scheme.
    //
    explicit KineticPlasmaDesc() noexcept = default;
    constexpr KineticPlasmaDesc(PlasmaDesc const &desc, unsigned Nc, ShapeOrder shape_order,
                                ParticleScheme scheme = full_f)
    : PlasmaDesc(desc), Nc{Nc}, shape_order{shape_order}, scheme{scheme}
    {
        if (this->Nc <= 0)
            throw std::invalid_argument{"Nc should be positive"};
    }

private:
    [[nodiscard]] friend constexpr auto serialize(KineticPlasmaDesc const &desc) noexcept
    {
        PlasmaDesc const &base = desc;
        return std::tuple_cat(serialize(base), std::make_tuple(desc.Nc, desc.scheme));
    }
};

/// Bi-Maxwellian plasma descriptor.
///
struct [[nodiscard]] BiMaxPlasmaDesc : public KineticPlasmaDesc {
    Real beta1; //!< The parallel component of plasma beta.
    Real T2_T1; //!< The ratio of the perpendicular to parallel temperatures.
    Real Vd;    //!< Equilibrium parallel drift speed.
    //
    constexpr BiMaxPlasmaDesc(KineticPlasmaDesc const &desc, Real beta1, Real T2_T1, Real Vd)
    : KineticPlasmaDesc(desc), beta1{beta1}, T2_T1{T2_T1}, Vd{Vd}
    {
        if (this->beta1 <= 0)
            throw std::invalid_argument{"beta1 should be positive"};
        if (this->T2_T1 <= 0)
            throw std::invalid_argument{"T2_T1 should be positive"};
    }
    constexpr BiMaxPlasmaDesc(KineticPlasmaDesc const &desc, Real beta1, Real T2_T1 = 1)
    : BiMaxPlasmaDesc(desc, beta1, T2_T1, {})
    {
    }

private:
    [[nodiscard]] friend constexpr auto serialize(BiMaxPlasmaDesc const &desc) noexcept
    {
        KineticPlasmaDesc const &base = desc;
        return std::tuple_cat(serialize(base), std::make_tuple(desc.beta1, desc.T2_T1, desc.Vd));
    }
};

/// Losscone distribution plasma descriptor.
///
/// The effective perpendicular temperature is 2*T2 = 1 + (1 - Δ)*β.
///
struct [[nodiscard]] LossconePlasmaDesc : public BiMaxPlasmaDesc {
    Real Delta; // Losscone VDF Δ parameter.
    Real beta;  // Losscone VDF β parameter.
    //
    constexpr LossconePlasmaDesc(BiMaxPlasmaDesc const &desc, Real Delta = 1, Real beta = 1)
    : BiMaxPlasmaDesc(desc), Delta{Delta}, beta{beta}
    {
        if (this->Delta < 0 || this->Delta > 1)
            throw std::invalid_argument{"Losscone.Delta should be in the range of [0, 1]"};
        if (this->beta <= 0)
            throw std::invalid_argument{"Losscone.beta should be positive"};
    }
    explicit constexpr LossconePlasmaDesc(KineticPlasmaDesc const &desc, Real beta1,
                                          Real vth_ratio /*ratio of θ2^2/θ1^2*/ = 1,
                                          std::pair<Real, Real> Db = {1, 1}, Real Vd = 0)
    : LossconePlasmaDesc({desc, beta1, (1 + (1 - Db.first) * Db.second) * vth_ratio, Vd}, Db.first,
                         Db.second)
    {
    }

private:
    [[nodiscard]] friend constexpr auto serialize(LossconePlasmaDesc const &desc) noexcept
    {
        BiMaxPlasmaDesc const &base = desc;
        return std::tuple_cat(serialize(base), std::make_tuple(desc.Delta, desc.beta));
    }
};
PIC1D_END_NAMESPACE

#endif /* PlasmaDesc_h */
