//
//  EField.h
//  hybrid_1d
//
//  Created by KYUNGGUK MIN on 1/15/19.
//  Copyright © 2019 Kyungguk Min & Kaijun Liu. All rights reserved.
//

#ifndef EField_h
#define EField_h

#include "../Geometry.h"
#include "../ParamSet.h"

HYBRID1D_BEGIN_NAMESPACE
class BField;
class Charge;
class Current;

class EField : public VectorGrid {
    VectorGrid Je;
    ScalarGrid Pe;

public:
    ParamSet const params;
    Geometry const geomtr;

public:
    explicit EField(ParamSet const &);

    void update(BField const &bfield, Charge const &charge, Current const &current) noexcept;

private:
    inline void _update_Pe(ScalarGrid &Pe, Charge const &rho) const noexcept;
    inline void _update_Je(VectorGrid &Je, Current const &Ji, BField const &B) const noexcept;
    inline void _update_E(EField &E, BField const &B, Charge const &rho) const noexcept;
};
HYBRID1D_END_NAMESPACE

#endif /* EField_h */
