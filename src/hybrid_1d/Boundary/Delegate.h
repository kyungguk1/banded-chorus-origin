//
//  Delegate.h
//  hybrid_1d
//
//  Created by KYUNGGUK MIN on 1/18/19.
//  Copyright © 2019 Kyungguk Min & Kaijun Liu. All rights reserved.
//

#ifndef Delegate_h
#define Delegate_h

#include "../Core/Domain.h"
#include "../InputWrapper.h"
#include "../Utility/GridQ.h"

HYBRID1D_BEGIN_NAMESPACE
class Delegate {
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
    using PartBucket = PartSpecies::bucket_type;
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
