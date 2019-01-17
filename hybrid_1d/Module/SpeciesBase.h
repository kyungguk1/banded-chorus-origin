//
//  SpeciesBase.h
//  hybrid_1d
//
//  Created by KYUNGGUK MIN on 1/17/19.
//  Copyright © 2019 kyungguk.com. All rights reserved.
//

#ifndef SpeciesBase_h
#define SpeciesBase_h

#include "../Utility/GridQ.h"
#include "../Utility/Scalar.h"
#include "../Utility/Vector.h"
#include "../Utility/Tensor.h"
#include "../Predefined.h"
#include "../Macros.h"
#include "../Inputs.h"

#include <tuple>
#include <vector>

HYBRID1D_BEGIN_NAMESPACE
class _Species {
public:
    // single particle description
    //
    struct Particle {
        Vector vel;
        Scalar pos;
    };

    // member variables
    //
    Real Nc; //!< number particles per cell
    Real Oc; //!< cyclotron frequency
    Real op; //!< plasma frequency
    std::vector<Particle> bucket; //!< particle container
private:
    std::tuple<GridQ<Scalar>, GridQ<Vector>, GridQ<Tensor>> _mom; //!< velocity moments at grid points

public:
    // accessors
    //
    Real charge_density_conversion_factor() const noexcept {
        return op*op*Input::O0/Oc;
    }
    Real current_density_conversion_factor() const noexcept {
        return charge_density_conversion_factor()/Input::c;
    }
    Real energy_density_conversion_factor() const noexcept {
        Real const tmp = Input::O0/Oc*op/Input::c;
        return tmp*tmp;
    }

    // access to i'th velocity moment
    //
    template <long i>
    auto moment() const noexcept -> std::tuple_element_t<i, decltype(_mom)> const &{
        return std::get<i>(_mom);
    }
    template <long i>
    auto moment()       noexcept -> std::tuple_element_t<i, decltype(_mom)>       &{
        return std::get<i>(_mom);
    }

protected:
    // constructor
    //
    explicit _Species(Real const Oc, Real const op, long const Nc);
};
HYBRID1D_END_NAMESPACE

#endif /* SpeciesBase_h */
