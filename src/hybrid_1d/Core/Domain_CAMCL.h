//
//  Domain_CAMCL.h
//  hybrid_1d
//
//  Created by KYUNGGUK MIN on 1/25/19.
//  Copyright © 2019 Kyungguk Min & Kaijun Liu. All rights reserved.
//

#ifndef Domain_CAMCL_h
#define Domain_CAMCL_h

#include "./Domain.h"

HYBRID1D_BEGIN_NAMESPACE
/// current-advanced cyclic leap-frog (CAM-CL) algorithm
///
class Domain_CAMCL : public Domain {
    // workspaces
    //
    BField  bfield_1;
    Current current_1;
    Charge  charge_1;
    Lambda  lambda;
    Gamma   gamma;

public:
    explicit Domain_CAMCL(ParamSet const &params, Delegate *delegate);

private:
    void advance_by(unsigned n_steps) override;
    void cycle(Domain const &domain);
    void subcycle(Domain const &, Charge const &charge, Current const &current, Real dt);
};
HYBRID1D_END_NAMESPACE

#endif /* Domain_CAMCL_h */
