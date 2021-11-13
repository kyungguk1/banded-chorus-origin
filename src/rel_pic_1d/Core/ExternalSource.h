/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include "Species.h"
#include <PIC/Badge.h>
#include <PIC/PlasmaDesc.h>

#include <HDF5Kit/HDF5Kit.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

PIC1D_BEGIN_NAMESPACE
class MasterDelegate;
class WorkerDelegate;

/// external current source
///
class ExternalSource : public Species {
    ExternalSourceBase      src_desc;
    std::vector<CurviCoord> src_pos;
    std::vector<Vector>     src_Jre;
    std::vector<Vector>     src_Jim;
    Real                    m_weighting_factor{ std::numeric_limits<Real>::quiet_NaN() };
    Real                    ramp_slope{ std::numeric_limits<Real>::quiet_NaN() };
    long                    m_cur_step{ std::numeric_limits<long>::quiet_NaN() };
    unsigned                number_of_source_points;

public:
    [[nodiscard]] PlasmaDesc const *operator->() const noexcept override { return &src_desc; }
    //
    [[nodiscard]] Real charge_density_conversion_factor() const noexcept override { return 1; }
    [[nodiscard]] Real current_density_conversion_factor() const noexcept override { return 1; }
    [[nodiscard]] Real energy_density_conversion_factor() const noexcept override { return 1; }

    // this is a hack to allow Master/WorkerDelegate to modify equilibrium_macro_weight
    [[nodiscard]] auto &weighting_factor(Badge<MasterDelegate>) &noexcept { return m_weighting_factor; }
    [[nodiscard]] auto &weighting_factor(Badge<WorkerDelegate>) &noexcept { return m_weighting_factor; }

    ExternalSource &operator=(ExternalSource const &) = delete;
    template <unsigned N>
    ExternalSource(ParamSet const &params, ExternalSourceDesc<N> const &src)
    : Species{ params }, src_desc{ src }, src_pos{ begin(src.pos), end(src.pos) }, src_Jre(N), src_Jim(N), number_of_source_points(N)
    {
        std::transform(begin(src.J0), end(src.J0), begin(src_Jre), [](auto const &cv) noexcept -> Vector {
            return { cv.x.real(), cv.y.real(), cv.z.real() };
        });
        std::transform(begin(src.J0), end(src.J0), begin(src_Jim), [](auto const &cv) noexcept -> Vector {
            return { cv.x.imag(), cv.y.imag(), cv.z.imag() };
        });

        // evenly divide up the source contribution among the distributed particle subdomain clones
        m_weighting_factor = Real{ 1 } / params.number_of_distributed_particle_subdomain_clones;

        // ramp slope
        constexpr auto eps = 1e-15;
        (ramp_slope = M_PI) /= src_desc.ease_in > eps ? src_desc.ease_in : 1.0;
    }
    ExternalSource() = default; // needed for empty std::array

    /// Reset the current step count
    /// \param step_count Simulation step count.
    void               set_cur_step(long step_count) noexcept { this->m_cur_step = step_count; }
    [[nodiscard]] long cur_step() const noexcept { return m_cur_step; }

    /// Update the moments at a given time
    /// \note After this call, the current step count is incremented.
    ///       Therefore, this should be called only once at each cycle.
    /// \param delta_t Time (delta from cur_step * params.dt) at which the moments should be calculated.
    void update(Real delta_t);

#ifndef DEBUG
private:
#endif
    void                 collect(VectorGrid &J, Real t) const;
    [[nodiscard]] Vector current(Vector const &J0re, Vector const &J0im, Real t) const noexcept;
    [[nodiscard]] Real   envelope(Real t) const noexcept;

    // attribute export facility
    //
    template <class Object>
    friend auto write_attr(Object &obj, ExternalSource const &sp) -> decltype(obj);
    friend auto operator<<(hdf5::Group &obj, ExternalSource const &sp) -> decltype(obj);
    friend auto operator<<(hdf5::Dataset &obj, ExternalSource const &sp) -> decltype(obj);
    friend auto operator<<(hdf5::Group &&obj, ExternalSource const &sp) -> decltype(obj)
    {
        return std::move(obj << sp);
    }
    friend auto operator<<(hdf5::Dataset &&obj, ExternalSource const &sp) -> decltype(obj)
    {
        return std::move(obj << sp);
    }
};
PIC1D_END_NAMESPACE
