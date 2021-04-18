//
//  VHistogramRecorder.h
//  hybrid_1d
//
//  Created by KYUNGGUK MIN on 4/15/20.
//  Copyright © 2020 Kyungguk Min & Kaijun Liu. All rights reserved.
//

#ifndef VHistogramRecorder_h
#define VHistogramRecorder_h

#include "./Recorder.h"

#include <fstream>
#include <string>

HYBRID1D_BEGIN_NAMESPACE
/// gyro-averaged velocity histogram recorder
///
/// particle samples over all domain are counted.
/// the histogram returned is normalized by the number of samples used to contruct the histogram
///
class VHistogramRecorder : public Recorder {
    std::ofstream os;

public:
    explicit VHistogramRecorder(unsigned rank, unsigned size);

private:
    std::string filepath(std::string const &wd, long step_count, unsigned sp_id) const;

    void record(Domain const &domain, long step_count) override;

    class Indexer;
    using vhist_t = std::map<std::pair<long, long>, std::pair<Real, Real>>;
    vhist_t histogram(PartSpecies const &sp, Indexer const &idxer) const;
    vhist_t histogram(Indexer const &idxer) const;
};
HYBRID1D_END_NAMESPACE

#endif /* VHistogramRecorder_h */
