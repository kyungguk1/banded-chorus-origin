/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include "Species.h"
#include <PIC/Badge.h>
#include <PIC/Particle.h>
#include <PIC/PlasmaDesc.h>
#include <PIC/VDFVariant.h>

#include <HDF5Kit/HDF5Kit.h>
#include <deque>
#include <sstream>
#include <vector>

HYBRID1D_BEGIN_NAMESPACE
class EField;
class BField;
class MasterDelegate;
class WorkerDelegate;

/// discrete simulation particle species
///
class PartSpecies : public Species {
    KineticPlasmaDesc desc;
    VDFVariant        vdf;
    Real              Nc;                         //!< number of particles per cell at the equator to be used for normalization
    Real              m_equilibrium_macro_weight; // weighting factor for delta-f equilibrium macro quantities
    ScalarGrid        equilibrium_mom0;
    VectorGrid        equilibrium_mom1;
    TensorGrid        equilibrium_mom2;

public:
    using bucket_type = std::deque<Particle>;
    bucket_type bucket; //!< particle container

    [[nodiscard]] KineticPlasmaDesc const *operator->() const noexcept override { return &desc; }

    // this is a hack to allow MasterDelegate to modify equilibrium_macro_weight
    [[nodiscard]] auto &equilibrium_macro_weight(Badge<MasterDelegate>) &noexcept { return m_equilibrium_macro_weight; }
    [[nodiscard]] auto &equilibrium_macro_weight(Badge<WorkerDelegate>) &noexcept { return m_equilibrium_macro_weight; }

    PartSpecies &operator=(PartSpecies const &) = default;
    PartSpecies(ParamSet const &params, KineticPlasmaDesc const &desc, VDFVariant vdf);
    PartSpecies() = default; // needed for empty std::array
    explicit PartSpecies(ParamSet const &params)
    : Species{ params } {} // needed for Domain_PC

    /// Load particles
    /// \note This should only be called by master thread.
    /// \param color This instructs which particles to keep.
    ///              Say, i'th particle is being loaded. It is kept if `color == i % divisor`.
    /// \param divisor The number of groups to which particles are divided.
    void populate(long color, long divisor);

    // load particles from a snapshot
    void load_ptls(std::vector<Particle> const &payload, bool append = false);

    // dump particles
    [[nodiscard]] std::vector<Particle> dump_ptls() const;

    void update_vel(BField const &bfield, EField const &efield, Real dt);
    void update_pos(Real dt, Real fraction_of_grid_size_allowed_to_travel);

    void collect_part(); // collect 0th and 1st moments
    void collect_all();  // collect all moments

private:
    void (PartSpecies::*m_update_velocity)(bucket_type &, VectorGrid const &, EField const &, BorisPush const &) const;
    void (PartSpecies::*m_collect_full_f)(ScalarGrid &, VectorGrid &) const;
    void (PartSpecies::*m_collect_delta_f)(ScalarGrid &, VectorGrid &, bucket_type &) const;

    [[nodiscard]] bool impl_update_pos(bucket_type &bucket, Real dt, Real travel_scale_factor) const;

    template <long Order>
    void impl_update_velocity(bucket_type &bucket, VectorGrid const &B, EField const &E, BorisPush const &boris) const;

    template <long Order>
    void impl_collect_full_f(ScalarGrid &, VectorGrid &nV) const; // weight is untouched
    template <long Order>
    void impl_collect_delta_f(ScalarGrid &, VectorGrid &nV, bucket_type &bucket) const; // weight is updated
    void impl_collect(ScalarGrid &n, VectorGrid &nV, TensorGrid &nvv) const;

    // attribute export facility
    //
    template <class Object>
    friend auto write_attr(Object &obj, PartSpecies const &sp) -> decltype(obj);
    friend auto operator<<(hdf5::Group &obj, PartSpecies const &sp) -> decltype(obj);
    friend auto operator<<(hdf5::Dataset &obj, PartSpecies const &sp) -> decltype(obj);
    friend auto operator<<(hdf5::Group &&obj, PartSpecies const &sp) -> decltype(obj)
    {
        return std::move(obj << sp);
    }
    friend auto operator<<(hdf5::Dataset &&obj, PartSpecies const &sp) -> decltype(obj)
    {
        return std::move(obj << sp);
    }
};
HYBRID1D_END_NAMESPACE
