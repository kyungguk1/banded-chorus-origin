//
//  BField.h
//  hybrid_1d
//
//  Created by KYUNGGUK MIN on 1/15/19.
//  Copyright © 2019 Kyungguk Min & Kaijun Liu. All rights reserved.
//

#ifndef BField_h
#define BField_h

#include "../Geometry.h"
#include "../ParamSet.h"

HYBRID1D_BEGIN_NAMESPACE
class EField;

class BField : public VectorGrid {
public:
    ParamSet const params;
    Geometry const geomtr;

public:
    explicit BField(ParamSet const &);
    BField &operator=(BField const &o) noexcept
    {
        this->GridQ::operator=(o);
        return *this;
    }
    using GridQ::swap;

    void update(EField const &efield, Real dt) noexcept;

private:
    static inline void _update(BField &B, EField const &E, Real cdtODx) noexcept;
};
HYBRID1D_END_NAMESPACE

#endif /* BField_h */
