/*
 * Copyright (c) 2020-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include "Recorder.h"

#include <map>
#include <string>

HYBRID1D_BEGIN_NAMESPACE
/// gyro-averaged velocity histogram recorder
///
/// particle samples over all domain are counted.
/// the histogram returned is normalized by the number of samples used to contruct the histogram
///
class VHistogramRecorder : public Recorder {
public:
    explicit VHistogramRecorder(parallel::mpi::Comm comm);

private:
    [[nodiscard]] std::string filepath(std::string const &wd, long step_count) const;

    class Indexer;
    using global_vhist_t = std::map<vhist_key_t, std::pair<Real, Real>>;
    using local_vhist_t  = std::map<vhist_key_t, vhist_val_t>;

    global_vhist_t histogram(PartSpecies const &sp, Indexer const &idxer) const;
    global_vhist_t histogram(Indexer const &idxer) const;

    void record(Domain const &domain, long step_count) override;
    void record_master(Domain const &domain, long step_count);
    void record_worker(Domain const &domain, long step_count);

    template <class Object>
    static decltype(auto) write_attr(Object &&obj, Domain const &domain, long step);
    static void           write_data(hdf5::Group &root, global_vhist_t vhist);
};
HYBRID1D_END_NAMESPACE
