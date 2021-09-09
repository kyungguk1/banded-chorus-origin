/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include "../ParamSet.h"
#include <PIC/BorisPush.h>

#include <HDF5Kit/HDF5Kit.h>
#include <tuple>

HYBRID1D_BEGIN_NAMESPACE
/// base class for ion species
///
class Species {
public:
    ParamSet const params;

protected:
    using MomTuple = std::tuple<ScalarGrid, VectorGrid, TensorGrid>;

private:
    MomTuple _mom{}; //!< velocity moments at grid points

    template <class T>
    using grid_t = Grid<T, ScalarGrid::size(), Pad>;

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
    template <long i>
    [[nodiscard]] auto const &moment() const noexcept
    {
        return std::get<i>(_mom);
    }
    template <long i>
    [[nodiscard]] auto &moment() noexcept { return std::get<i>(_mom); }
    //
    template <class T>
    [[nodiscard]] auto const &moment() const noexcept
    {
        return std::get<grid_t<T>>(_mom);
    }
    template <class T>
    [[nodiscard]] auto &moment() noexcept { return std::get<grid_t<T>>(_mom); }
    //
    [[nodiscard]] MomTuple const &moments() const noexcept { return _mom; }
    [[nodiscard]] MomTuple       &moments() noexcept { return _mom; }

protected:
    virtual ~Species() = default;
    Species(ParamSet const & = {});
    Species &operator=(Species const &) noexcept;
    Species &operator=(Species &&) noexcept;

    // attribute export facility
    //
    friend auto operator<<(hdf5::Group &obj, Species const &sp) -> decltype(obj);
    friend auto operator<<(hdf5::Dataset &obj, Species const &sp) -> decltype(obj);
    friend auto operator<<(hdf5::Group &&obj, Species const &sp) -> decltype(obj)
    {
        return std::move(obj << sp);
    }
    friend auto operator<<(hdf5::Dataset &&obj, Species const &sp) -> decltype(obj)
    {
        return std::move(obj << sp);
    }
};
HYBRID1D_END_NAMESPACE
