/*
 * Copyright (c) 2020-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include "Recorder.h"

#include <map>
#include <string>

PIC1D_BEGIN_NAMESPACE
/// gyro-averaged velocity histogram recorder
///
/// particle samples over all domain are counted.
/// the histogram returned is normalized by the number of samples used to contruct the histogram
///
class VHistogramRecorder : public Recorder {
    class Indexer;
    using index_pair_t   = std::pair<long, long>;
    using local_vhist_t  = std::map<index_pair_t, std::pair<unsigned long, Real>>;
    using global_vhist_t = std::map<index_pair_t, std::pair<Real, Real>>;

    parallel::Communicator<unsigned long, local_vhist_t::value_type> world_comm;

public:
    VHistogramRecorder(parallel::mpi::Comm subdomain_comm, parallel::mpi::Comm const &world_comm);

private:
    [[nodiscard]] std::string filepath(std::string const &wd, long step_count) const;

    void record(Domain const &domain, long step_count) override;
    void record_master(Domain const &domain, long step_count);
    void record_worker(Domain const &domain, long step_count);

    auto               histogram(PartSpecies const &sp, Indexer const &idxer) const -> global_vhist_t;
    [[nodiscard]] auto local_counting(PartSpecies const &sp, Indexer const &idxer) const -> local_vhist_t;
    [[nodiscard]] auto global_counting(unsigned long local_count, local_vhist_t local_vhist) const
        -> std::pair<unsigned long /*total count*/, local_vhist_t>;

    template <class Object>
    static decltype(auto) write_attr(Object &&obj, Domain const &domain, long step);
    static void           write_data(hdf5::Group &root, global_vhist_t vhist);
};
PIC1D_END_NAMESPACE
