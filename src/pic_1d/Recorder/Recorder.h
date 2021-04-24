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

#ifndef Recorder_h
#define Recorder_h

#include "../Core/Domain.h"
#include "../TypeMaps.h"

#include <ParallelKit/ParallelKit.h>

PIC1D_BEGIN_NAMESPACE
class Recorder {
public:
    virtual ~Recorder()                                        = default;
    virtual void record(Domain const &domain, long step_count) = 0;

    using vhist_key_t     = std::pair<long, long>; // {v1, v2} indices
    using vhist_val_t     = std::pair<long, Real>; // {full-f, delta-f} vhist
    using vhist_payload_t = std::pair<vhist_key_t const, vhist_val_t>;
    using interprocess_comm_t
        = parallel::Communicator<Scalar, Vector, Tensor, Particle, vhist_payload_t,
                                 unsigned long /*local particle count*/>;
    using rank_t = parallel::mpi::Rank;

public:
    long const recording_frequency;

protected:
    interprocess_comm_t const comm;

    static constexpr parallel::mpi::Tag tag{875};
    static constexpr char               null_dev[] = "/dev/null";
    static constexpr rank_t             master{0};
    [[nodiscard]] bool                  is_master() const { return master == comm->rank(); }

protected:
    explicit Recorder(unsigned recording_frequency, parallel::mpi::Comm comm);
};
PIC1D_END_NAMESPACE

#endif /* Recorder_h */
