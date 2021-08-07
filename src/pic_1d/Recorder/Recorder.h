/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#ifndef Recorder_h
#define Recorder_h

#include "../Core/Domain.h"
#include "../Utility/TypeMaps.h"

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

    static constexpr parallel::mpi::Tag tag{ 875 };
    static constexpr char               null_dev[] = "/dev/null";
    static constexpr rank_t             master{ 0 };
    [[nodiscard]] bool                  is_master() const { return master == comm->rank(); }

protected:
    explicit Recorder(unsigned recording_frequency, parallel::mpi::Comm comm);
};
PIC1D_END_NAMESPACE

#endif /* Recorder_h */
