//
//  MomentRecorder.h
//  pic_1d
//
//  Created by KYUNGGUK MIN on 1/29/19.
//  Copyright © 2019 Kyungguk Min & Kaijun Liu. All rights reserved.
//

#ifndef MomentRecorder_h
#define MomentRecorder_h

#include "./Recorder.h"

#include <fstream>
#include <string>

PIC1D_BEGIN_NAMESPACE
/// ion moment recorder
/// field-aligned components are recorded;
/// suffix 1, 2, and 3 means three field-aligned components:
///     1 : parallel, 2 : perpendicular, and 3 : out-of-plane
///
class MomentRecorder : public Recorder {
    std::ofstream os;

public:
    explicit MomentRecorder(unsigned rank, unsigned size);

private:
    std::string filepath(std::string const &wd, long step_count) const;

    void record(Domain const &domain, long step_count) override;
};
PIC1D_END_NAMESPACE

#endif /* MomentRecorder_h */
