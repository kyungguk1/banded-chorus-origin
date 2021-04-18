//
//  ParticleRecorder.h
//  pic_1d
//
//  Created by KYUNGGUK MIN on 1/29/19.
//  Copyright © 2019 Kyungguk Min & Kaijun Liu. All rights reserved.
//

#ifndef ParticleRecorder_h
#define ParticleRecorder_h

#include "./Recorder.h"

#include <fstream>
#include <random>
#include <string>

PIC1D_BEGIN_NAMESPACE
/// marker particle recorder
/// field-aligned components are recorded;
/// suffix 1, 2, and 3 means three field-aligned components:
///     1 : parallel, 2 : perpendicular, and 3 : out-of-plane
///
class ParticleRecorder : public Recorder {
    std::mt19937  urbg;
    std::ofstream os;

public:
    explicit ParticleRecorder(unsigned rank, unsigned size);

private:
    std::string filepath(std::string const &wd, long step_count, unsigned sp_id) const;

    void record(Domain const &domain, long step_count) override;
    void record(PartSpecies const &sp, unsigned max_count);
};
PIC1D_END_NAMESPACE

#endif /* ParticleRecorder_h */
