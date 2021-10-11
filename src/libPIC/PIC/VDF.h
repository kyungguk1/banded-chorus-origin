/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC/Config.h>
#include <PIC/CurviCoord.h>
#include <PIC/Geometry.h>
#include <PIC/Particle.h>
#include <PIC/PlasmaDesc.h>
#include <PIC/Predefined.h>
#include <PIC/Range.h>
#include <PIC/Scalar.h>
#include <PIC/Tensor.h>
#include <PIC/Vector.h>

#include <vector>

LIBPIC_BEGIN_NAMESPACE
/// Base class for velocity distribution function
///
template <class Concrete>
class VDF {
    using Self = Concrete;

    [[nodiscard]] constexpr decltype(auto) self() const noexcept { return static_cast<Self const *>(this); }
    [[nodiscard]] constexpr decltype(auto) self() noexcept { return static_cast<Self *>(this); }

protected:
    Geometry geomtr;
    Range    domain_extent; // particle boundary extent

    VDF(Geometry const &geo, Range const &domain_extent)
    noexcept
    : geomtr{ geo }, domain_extent{ domain_extent }
    {
    }

public:
    /// Plasma description associated with *this
    ///
    [[nodiscard]] decltype(auto) plasma_desc() const noexcept { return self()->impl_plasma_desc(); }

    /// Sample a single particle following the marker particle distribution, g0.
    /// \note Concrete subclass should provide impl_emit with the same signature.
    ///
    [[nodiscard]] Particle emit() const { return self()->impl_emit(); }

    /// Sample N particles following the marker particle distribution, g0.
    /// \note Concrete subclass should provide impl_emit with the same signature.
    [[nodiscard]] std::vector<Particle> emit(unsigned long n) const { return self()->impl_emit(n); }

    /// Zero velocity moment at the given position, \<1\>_0(x).
    /// \note Concrete subclass should provide impl_n with the same signature.
    ///
    [[nodiscard]] Scalar n0(CurviCoord const &pos) const { return self()->impl_n(pos); }

    /// First velocity moment at the given position, \<v\>_0(x).
    /// \note Concrete subclass should provide impl_nV with the same signature.
    ///
    [[nodiscard]] Vector nV0(CurviCoord const &pos) const { return self()->impl_nV(pos); }

    /// Second velocity moment at the given position, \<vv\>_0(x).
    /// \note Concrete subclass should provide impl_nvv with the same signature.
    ///
    [[nodiscard]] Tensor nvv0(CurviCoord const &pos) const { return self()->impl_nvv(pos); }

    /// Calculate delta-f weighting factor
    /// \details The weight is given by
    ///
    /// f(0, x(0), v(0))/g(0, x(0), v(0)) - f_0(x(t), v(t))/g(0, x(0), v(0))
    ///
    /// where g is the marker particle distribution.
    /// \note Concrete subclass should provide impl_weight with the same signature.
    ///
    [[nodiscard]] Real weight(Particle const &ptl) const { return self()->impl_weight(ptl); }

    /// Ratio of the number of particles at the reference cell to the total number of particles
    /// \note Concrete subclass should provide a member variable, m_Nrefcell_div_Ntotal, containing this quantity.
    [[nodiscard]] Real Nrefcell_div_Ntotal() const { return self()->m_Nrefcell_div_Ntotal; }
};
LIBPIC_END_NAMESPACE
