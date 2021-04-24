/*
 * Copyright (c) 2020-2021, Kyungguk Min
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

#ifndef VHistogramRecorder_h
#define VHistogramRecorder_h

#include "./Recorder.h"

#include <fstream>
#include <map>
#include <string>

PIC1D_BEGIN_NAMESPACE
/// gyro-averaged velocity histogram recorder
///
/// particle samples over all domain are counted.
/// the histogram returned is normalized by the number of samples used to contruct the histogram
///
class VHistogramRecorder : public Recorder {
    std::ofstream os;

public:
    explicit VHistogramRecorder(parallel::mpi::Comm comm);

private:
    std::string filepath(std::string const &wd, long step_count, unsigned sp_id) const;

    void record(Domain const &domain, long step_count) override;

    class Indexer;
    using global_vhist_t = std::map<vhist_key_t, std::pair<Real, Real>>;
    using local_vhist_t  = std::map<vhist_key_t, vhist_val_t>;

    global_vhist_t histogram(PartSpecies const &sp, Indexer const &idxer) const;
    global_vhist_t histogram(Indexer const &idxer) const;
};
PIC1D_END_NAMESPACE

#endif /* VHistogramRecorder_h */
