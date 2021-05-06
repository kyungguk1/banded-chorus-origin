/*
 * Copyright (c) 2019-2021, Kyungguk Min
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * - Redistributions of source code must retain the above copyright notice, this
 *  list of conditions and the following disclaimer.
 *
 * - Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef PartSpecies_h
#define PartSpecies_h

#include "../VDF/VDF.h"
#include "./Species.h"

#include <HDF5Kit/HDF5Kit.h>
#include <deque>
#include <memory>
#include <sstream>
#include <vector>

PIC1D_BEGIN_NAMESPACE
class EField;
class BField;

/// discrete simulation particle species
///
class PartSpecies : public Species {
    KineticPlasmaDesc    desc;
    std::unique_ptr<VDF> vdf;

public:
    using bucket_type = std::deque<Particle>;
    bucket_type bucket; //!< particle container
    Real Nc; //!< number of particles per cell to be used to normalization; don't modify this if you
             //!< don't know what you are doing

public:
    [[nodiscard]] KineticPlasmaDesc const *operator->() const noexcept override { return &desc; }

    PartSpecies &operator=(PartSpecies &&) = delete;
    PartSpecies()                          = default; // needed for empty std::array
    PartSpecies(ParamSet const &params, KineticPlasmaDesc const &desc,
                std::unique_ptr<VDF> vdf); // leaves bucket empty

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
    void (*_update_velocity)(bucket_type &, VectorGrid const &, EField const &, BorisPush const);
    void (PartSpecies::*_collect_full_f)(VectorGrid &) const;
    void (PartSpecies::*_collect_delta_f)(VectorGrid &, bucket_type &) const;

private:
    [[nodiscard]] static bool _update_x(bucket_type &bucket, Real dtODx, Real travel_scale_factor);

    template <long Order>
    static void _update_velocity_(bucket_type &bucket, VectorGrid const &B, EField const &E,
                                  BorisPush pusher);

    template <long Order> void _collect_full_f_(VectorGrid &nV) const; // weight is untouched
    template <long Order>
    void _collect_delta_f_(VectorGrid &nV, bucket_type &bucket) const; // weight is updated
    void _collect(ScalarGrid &n, VectorGrid &nV, TensorGrid &nvv) const;

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
                          PartSpecies::bucket_type const &   bucket)
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

#endif /* PartSpecies_h */
