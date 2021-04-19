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
