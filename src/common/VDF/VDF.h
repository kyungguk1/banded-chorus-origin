/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <Geometry.h>
#include <PlasmaDesc.h>
#include <Predefined.h>
#include <Utility/Range.h>
#include <Utility/Scalar.h>
#include <Utility/Tensor.h>
#include <Utility/Vector.h>
#include <VDF/Particle.h>
#include <common-config.h>

#include <type_traits>
#include <vector>

COMMON_BEGIN_NAMESPACE
/// Base class for velocity distribution function
///
template <class Concrete> class VDF {
    using Self = Concrete;

    [[nodiscard]] constexpr auto self() const noexcept { return static_cast<Self const *>(this); }
    [[nodiscard]] constexpr auto self() noexcept { return static_cast<Self *>(this); }

protected:
    Geometry geomtr;
    Range    domain_extent;

    VDF(Geometry const &geo, Range const &domain_extent)
    noexcept : geomtr{ geo }, domain_extent{ domain_extent }
    {
    }

public:
    /// Convenience ctor of concrete VDF object
    ///
    template <class... Args> [[nodiscard]] static auto make(Args &&...args)
    {
        return Concrete{ std::forward<Args>(args)... };
    }

    /// Sample a single particle following the marker particle distribution, g0.
    /// \note Concrete subclass should provide impl_emit with the same signature.
    ///
    [[nodiscard]] Particle emit() const { return self()->impl_emit(); }

    /// Sample N particles following the marker particle distribution, g0.
    [[nodiscard]] auto emit(unsigned n) const
    {
        std::vector<Particle> ptls(n);
        for (auto &ptl : ptls)
            ptl = this->emit();
        return ptls;
    }

    /// Zero velocity moment at the given position, \<1\>_0(x).
    /// \note Concrete subclass should provide impl_n0 with the same signature.
    ///
    [[nodiscard]] Scalar n0(Real pos_x) const { return self()->impl_n0(pos_x); }

    /// First velocity moment at the given position, \<v\>_0(x).
    /// \note Concrete subclass should provide impl_nV0 with the same signature.
    ///
    [[nodiscard]] Vector nV0(Real pos_x) const { return self()->impl_nV0(pos_x); }

    /// Second velocity moment at the given position, \<vv\>_0(x).
    /// \note Concrete subclass should provide impl_nvv0 with the same signature.
    ///
    [[nodiscard]] Tensor nvv0(Real pos_x) const { return self()->impl_nvv0(pos_x); }

    /// Calculate the change of PSD
    /// \details Given a particle at some time instant, t,
    /// it calculate the change of PSD from the initial one:
    ///         1 - f_0(x(t), v(t))/f(0, x(0), v(0))
    ///
    /// \note Concrete subclass should provide impl_delta_f with the same signature.
    ///
    [[nodiscard]] Real delta_f(Particle const &ptl) const { return self()->impl_delta_f(ptl); }

    /// Calculate delta-f weighting factor
    /// \details The weight is given by
    ///
    /// f(0, x(0), v(0))/g(0, x(0), v(0)) - f_0(x(t), v(t))/g(0, x(0), v(0))
    ///
    /// where g is the marker particle distribution.
    ///
    [[nodiscard]] Real weight(Particle const &ptl) const { return ptl.psd.fOg * delta_f(ptl); }
};
COMMON_END_NAMESPACE
