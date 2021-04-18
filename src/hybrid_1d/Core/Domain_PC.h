//
//  Domain_PC.h
//  hybrid_1d
//
//  Created by KYUNGGUK MIN on 1/18/19.
//  Copyright © 2019 Kyungguk Min & Kaijun Liu. All rights reserved.
//

#ifndef Domain_PC_h
#define Domain_PC_h

#include "./Domain.h"

HYBRID1D_BEGIN_NAMESPACE
/// predictor-corrector algorithm
///
class Domain_PC : public Domain {
    // workspaces
    //
    BField      bfield_1;
    EField      efield_1;
    PartSpecies part_predict;
    ColdSpecies cold_predict;

public:
    explicit Domain_PC(ParamSet const &params, Delegate *delegate);

private:
    void advance_by(unsigned n_steps) override;
    void cycle(Domain const &domain);
    void predictor_step(Domain const &domain);
    void corrector_step(Domain const &domain);
};
HYBRID1D_END_NAMESPACE

#endif /* Domain_PC_h */
