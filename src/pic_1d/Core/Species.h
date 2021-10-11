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

PIC1D_BEGIN_NAMESPACE
/// base class for ion/electron species
///
class Species {
    template <class T>
    using grid_t   = Grid<T, ScalarGrid::size(), Pad>;
    using MomTuple = std::tuple<ScalarGrid, VectorGrid, TensorGrid>;

public:
    ParamSet const params;
    Geometry const geomtr;

private:
    MomTuple m_mom{}; //!< velocity moments at grid points

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
        return std::get<i>(m_mom);
    }
    template <long i>
    [[nodiscard]] auto &moment() noexcept { return std::get<i>(m_mom); }
    //
    template <class T>
    [[nodiscard]] auto const &moment() const noexcept
    {
        return std::get<grid_t<T>>(m_mom);
    }
    template <class T>
    [[nodiscard]] auto &moment() noexcept { return std::get<grid_t<T>>(m_mom); }
    //
    [[nodiscard]] MomTuple const &moments() const noexcept { return m_mom; }
    [[nodiscard]] MomTuple       &moments() noexcept { return m_mom; }

protected:
    virtual ~Species() = default;
    explicit Species(ParamSet const & = {});
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
PIC1D_END_NAMESPACE
