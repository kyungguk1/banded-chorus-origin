//
// Copyright (c) 2019, Kyungguk Min
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// - Redistributions of source code must retain the above copyright notice, this
//  list of conditions and the following disclaimer.
//
// - Redistributions in binary form must reproduce the above copyright notice,
//  this list of conditions and the following disclaimer in the documentation
//  and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#ifndef Recorder_h
#define Recorder_h

#include "../Core/Domain.h"
#include "../Utility/MessageDispatch.h"

HYBRID1D_BEGIN_NAMESPACE
class Recorder {
public:
    virtual ~Recorder()                                        = default;
    virtual void record(Domain const &domain, long step_count) = 0;

    using PartBucket = PartSpecies::bucket_type;
    using vhist_payload_t
        = std::map<std::pair<long, long>, std::pair<long, Real>>; // {full-f, delta-f} vhist
    using message_dispatch_t
        = MessageDispatch<Scalar, Vector, Tensor, PartBucket,
                          std::pair<Vector const *, Vector const *>,
                          std::pair<PartSpecies const *, ColdSpecies const *>,
                          std::pair<unsigned long /*local particle count*/, vhist_payload_t>>;
    using interthread_comm_t = message_dispatch_t::Communicator;

    long const                recording_frequency;
    std::vector<unsigned>     all_ranks;
    std::vector<unsigned>     all_but_master;
    static message_dispatch_t dispatch;
    interthread_comm_t const  comm;
    unsigned const            size;
    static constexpr char     null_dev[] = "/dev/null";
    static constexpr unsigned master     = 0;
    [[nodiscard]] bool        is_master() const noexcept { return master == comm.rank(); }

protected:
    explicit Recorder(unsigned recording_frequency, unsigned rank, unsigned size);

    template <class T, class Op> T reduce(T x, Op op)
    {
        if (is_master()) {
            // reduce; skip collecting master's value, cuz it is used as initial value
            x = comm.reduce<T>(all_but_master, x, op);
            // broadcast result
            auto tks = comm.bcast(x, all_but_master);
            return x;
            // std::for_each(std::make_move_iterator(begin(tks)), std::make_move_iterator(end(tks)),
            // std::mem_fn(&ticket_t::wait));
        } else {
            comm.send(x, master)
                .wait(); // wait can help to break the contention at the recv that follows
            return comm.recv<T>(master);
        }
    }
};
HYBRID1D_END_NAMESPACE

#endif /* Recorder_h */
