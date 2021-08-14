/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include "Species.h"
#include <PIC/PlasmaDesc.h>
#include <PIC/VDFVariant.h>

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
    KineticPlasmaDesc desc;
    VDFVariant        vdf;

public:
    [[nodiscard]] KineticPlasmaDesc const *operator->() const noexcept override { return &desc; }

    using bucket_type = std::deque<Particle>;
    bucket_type bucket; //!< particle container
    Real Nc; //!< number of particles per cell to be used to normalization; don't modify this if you
             //!< don't know what you are doing

    PartSpecies &operator=(PartSpecies &&) = delete;
    PartSpecies(ParamSet const &params, KineticPlasmaDesc const &desc, VDFVariant const &vdf);
    PartSpecies() = default; // needed for empty std::array

    // load particles using VDF; should only be called by master thread
    void populate();

    // load particles from a snapshot; particles' coordinates are
    // expected to be relative to the whole domain
    void load_ptls(std::vector<Particle> const &payload, bool append = false);

    // dump particles whose coordinates are relative to the whole domain
    [[nodiscard]] std::vector<Particle> dump_ptls() const;

    void update_vel(BField const &bfield, EField const &efield, Real dt);
    void update_pos(Real dt, Real fraction_of_grid_size_allowed_to_travel);

    void collect_part(); // collect 1st moment
    void collect_all();  // collect all moments

private:
    void (*m_update_velocity)(bucket_type &, VectorGrid const &, EField const &, BorisPush const &);
    void (PartSpecies::*m_collect_full_f)(VectorGrid &) const;
    void (PartSpecies::*m_collect_delta_f)(VectorGrid &, bucket_type &) const;

private:
    [[nodiscard]] static bool impl_update_x(bucket_type &bucket, Real dtODx,
                                            Real travel_scale_factor);

    template <long Order>
    static void impl_update_velocity(bucket_type &bucket, VectorGrid const &B, EField const &E,
                                     BorisPush const &boris);

    template <long Order> void impl_collect_full_f(VectorGrid &nV) const; // weight is untouched
    template <long Order>
    void impl_collect_delta_f(VectorGrid &nV, bucket_type &bucket) const; // weight is updated
    void impl_collect(ScalarGrid &n, VectorGrid &nV, TensorGrid &nvv) const;

    // attribute export facility
    //
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

// MARK:- pretty print for particle container
//
template <class CharT, class Traits>
decltype(auto) operator<<(std::basic_ostream<CharT, Traits> &os,
                          PartSpecies::bucket_type const    &bucket)
{
    std::basic_ostringstream<CharT, Traits> ss;
    {
        ss.flags(os.flags());
        ss.imbue(os.getloc());
        ss.precision(os.precision());
        //
        auto it = bucket.cbegin(), end = bucket.cend();
        ss << '{';
        if (it != end) { // check if bucket is empty
            ss << *it++;
        }
        while (it != end) {
            ss << ", " << *it++;
        }
        ss << '}';
    }
    return os << ss.str();
}
PIC1D_END_NAMESPACE
