/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#ifndef Delegate_h
#define Delegate_h

#include "../Core/Domain.h"
#include "../InputWrapper.h"
#include "../Utility/GridQ.h"
#include "../Utility/Particle.h"

#include <vector>

HYBRID1D_BEGIN_NAMESPACE
class Delegate {
protected:
    using PartBucket = std::vector<Particle>;
    mutable struct bucket_pair { // be careful not to access it from multiple threads
        PartBucket L{};
        PartBucket R{};

        [[nodiscard]] decltype(auto) cleared()
        {
            L.clear();
            R.clear();
            return *this;
        }
    } buckets{}; // be sure to clear the contents before use

public:
    Delegate &operator=(Delegate const &) = delete;
    Delegate(Delegate const &)            = delete;
    virtual ~Delegate()                   = default;
    explicit Delegate() noexcept          = default;

    // all virtual's called by Domain are const qualified to remind that changing the state of this
    // during concurrent calls likely cause the race condition and other side effects

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
    virtual void pass(Domain const &, PartBucket &L_bucket, PartBucket &R_bucket) const;
    virtual void pass(Domain const &, PartSpecies &) const;
    virtual void pass(Domain const &, ColdSpecies &) const   = 0;
    virtual void pass(Domain const &, BField &) const        = 0;
    virtual void pass(Domain const &, EField &) const        = 0;
    virtual void pass(Domain const &, Charge &) const        = 0;
    virtual void pass(Domain const &, Current &) const       = 0;
    virtual void gather(Domain const &, Charge &) const      = 0;
    virtual void gather(Domain const &, Current &) const     = 0;
    virtual void gather(Domain const &, PartSpecies &) const = 0;

private: // helpers
    template <class T, long N> static void pass(GridQ<T, N> &);
    template <class T, long N> static void gather(GridQ<T, N> &);
};
HYBRID1D_END_NAMESPACE

#endif /* Delegate_h */
