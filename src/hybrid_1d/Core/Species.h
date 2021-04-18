/*
 * Copyright (c) 2019, Kyungguk Min
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * - Redistributions of source code must retain the above copyright notice, this
 *  list of conditions and the following disclaimer.
 *
 * - Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef Species_h
#define Species_h

#include "../Geometry.h"
#include "../ParamSet.h"
#include "../Utility/BorisPush.h"
#include "../Utility/Particle.h"

#include <tuple>

HYBRID1D_BEGIN_NAMESPACE
/// base class for ion species
///
class Species {
public:
    ParamSet const params;
    Geometry const geomtr;

protected:
    using MomTuple = std::tuple<ScalarGrid, VectorGrid, TensorGrid>;

private:
    MomTuple _mom{}; //!< velocity moments at grid points

    template <class T> using grid_t = GridQ<T, ScalarGrid::size()>;

public:
    // accessors
    //
    [[nodiscard]] virtual PlasmaDesc const *operator->() const noexcept = 0;

    [[nodiscard]] Real charge_density_conversion_factor() const noexcept
    {
        return ((*this)->op * (*this)->op) * params.O0 / (*this)->Oc;
    }
    [[nodiscard]] Real current_density_conversion_factor() const noexcept
    {
        return charge_density_conversion_factor() / params.c;
    }
    [[nodiscard]] Real energy_density_conversion_factor() const noexcept
    {
        Real const tmp = params.O0 / (*this)->Oc * (*this)->op / params.c;
        return tmp * tmp;
    }

    // access to i'th velocity moment
    //
    template <long i> [[nodiscard]] auto const &moment() const noexcept
    {
        return std::get<i>(_mom);
    }
    template <long i> [[nodiscard]] auto &moment() noexcept { return std::get<i>(_mom); }
    //
    template <class T> [[nodiscard]] auto const &moment() const noexcept
    {
        return std::get<grid_t<T>>(_mom);
    }
    template <class T> [[nodiscard]] auto &moment() noexcept { return std::get<grid_t<T>>(_mom); }
    //
    [[nodiscard]] MomTuple const &moments() const noexcept { return _mom; }
    [[nodiscard]] MomTuple &      moments() noexcept { return _mom; }

protected:
    virtual ~Species() = default;
    explicit Species(ParamSet const & = {});
    Species &operator=(Species const &) noexcept;
    Species &operator=(Species &&) noexcept;
};
HYBRID1D_END_NAMESPACE

#endif /* Species_h */
