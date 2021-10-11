/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include "Species.h"
#include <PIC/PlasmaDesc.h>
#include <PIC/RelativisticParticle.h>
#include <PIC/RelativisticVDFVariant.h>

#include <HDF5Kit/HDF5Kit.h>
#include <deque>
#include <sstream>
#include <vector>

PIC1D_BEGIN_NAMESPACE
class EField;
class BField;

/// discrete simulation particle species
///
class PartSpecies : public Species {
    KineticPlasmaDesc      desc;
    RelativisticVDFVariant vdf;
    Real                   Nc;                       //!< number of particles per cell at the equator to be used for normalization
    Real                   equilibrium_macro_weight; // weighting factor for delta-f equilibrium macro quantities

    using Particle = RelativisticParticle;

public:
    using bucket_type = std::deque<Particle>;
    bucket_type bucket; //!< particle container

    [[nodiscard]] KineticPlasmaDesc const *operator->() const noexcept override { return &desc; }

    PartSpecies &operator=(PartSpecies &&) = delete; // this should not be default-ed
    PartSpecies(ParamSet const &params, KineticPlasmaDesc const &desc, RelativisticVDFVariant vdf);
    PartSpecies() = default; // needed for empty std::array

    // load particles using VDF; should only be called by master thread
    void populate();

    // load particles from a snapshot
    void load_ptls(std::vector<Particle> const &payload, bool append = false);

    // dump particles
    [[nodiscard]] std::vector<Particle> dump_ptls() const;

    void update_vel(BField const &bfield, EField const &efield, Real dt);
    void update_pos(Real dt, Real fraction_of_grid_size_allowed_to_travel);

    void collect_part(); // collect 1st moment; usual velocity moment
    void collect_all();  // collect all moments; moments of γv, i.e., relativistic momentum

private:
    void (PartSpecies::*m_update_velocity)(bucket_type &, VectorGrid const &, EField const &, BorisPush const &) const;
    void (PartSpecies::*m_collect_full_f)(VectorGrid &) const;
    void (PartSpecies::*m_collect_delta_f)(VectorGrid &, bucket_type &) const;

    [[nodiscard]] bool impl_update_pos(bucket_type &bucket, Real dt, Real travel_scale_factor) const;

    template <long Order>
    void impl_update_velocity(bucket_type &bucket, VectorGrid const &B, EField const &E, BorisPush const &boris) const;

    template <long Order>
    void impl_collect_full_f(VectorGrid &nV) const; // weight is untouched
    template <long Order>
    void impl_collect_delta_f(VectorGrid &nV, bucket_type &bucket) const; // weight is updated
    void impl_collect(ScalarGrid &n, VectorGrid &nV, FourTensorGrid &nuv) const;

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
PIC1D_END_NAMESPACE
