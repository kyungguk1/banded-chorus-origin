/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include "../Core/Domain.h"
#include "../ParamSet.h"
#include <PIC/Grid.h>
#include <PIC/RelativisticParticle.h>

#include <vector>

PIC1D_BEGIN_NAMESPACE
class Delegate {
protected:
    using Particle   = RelativisticParticle;
    using PartBucket = std::vector<Particle>;

    struct bucket_pair_t { // be careful not to access it from multiple threads
        PartBucket L{};
        PartBucket R{};

        [[nodiscard]] decltype(auto) cleared()
        {
            L.clear();
            R.clear();
            return (*this);
        }
    };
    mutable bucket_pair_t buckets{}; // be sure to clear the contents before use

public:
    Delegate &operator=(Delegate const &) = delete;
    Delegate(Delegate const &)            = delete;
    virtual ~Delegate()                   = default;
    explicit Delegate() noexcept          = default;

    // called once after initialization but right before entering loop
    //
    virtual void once(Domain &) const = 0;

    // called before and after every cycle of update
    //
    virtual void prologue(Domain const &, long inner_step_count) const = 0;
    virtual void epilogue(Domain const &, long inner_step_count) const = 0;

    // boundary value communication
    //
    virtual void partition(PartSpecies &, PartBucket &L_bucket, PartBucket &R_bucket) const;
    virtual void boundary_pass(Domain const &, PartBucket &L_bucket, PartBucket &R_bucket) const;
    virtual void boundary_pass(Domain const &, PartSpecies &) const; // be aware of mutation of particle bucket occurring
    virtual void boundary_pass(Domain const &, ColdSpecies &) const = 0;
    virtual void boundary_pass(Domain const &, BField &) const      = 0;
    virtual void boundary_pass(Domain const &, EField &) const      = 0;
    virtual void boundary_pass(Domain const &, Current &) const     = 0;
    virtual void boundary_gather(Domain const &, Current &) const   = 0;
    virtual void boundary_gather(Domain const &, Species &) const   = 0;

private: // helpers
    template <class T, long N>
    static void pass(Grid<T, N, Pad> &);
    template <class T, long N>
    static void gather(Grid<T, N, Pad> &);
    static void periodic_particle_pass(ParamSet const &, PartBucket &L_bucket, PartBucket &R_bucket);
    static void reflecting_particle_pass(ParamSet const &, PartBucket &L_bucket, PartBucket &R_bucket);

    [[nodiscard]] static Vector randomize_gyrophase(Vector const &);
};
PIC1D_END_NAMESPACE
