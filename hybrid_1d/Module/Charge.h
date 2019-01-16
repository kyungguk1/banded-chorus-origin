//
//  Charge.h
//  hybrid_1d
//
//  Created by KYUNGGUK MIN on 1/15/19.
//  Copyright © 2019 kyungguk.com. All rights reserved.
//

#ifndef Charge_h
#define Charge_h

#include "../Inputs.h"
#include "../Utility/GridQ.h"
#include "../Utility/Scalar.h"

HYBRID1D_BEGIN_NAMESPACE
class BField;
class EField;
class Species;

class Charge : protected GridQ<Scalar> {
public:
};
HYBRID1D_END_NAMESPACE

#endif /* Charge_h */