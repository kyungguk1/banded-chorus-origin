/*
 * Copyright (c) 2019-2022, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include "../ParamSet.h"

#include <HDF5Kit/HDF5Kit.h>

PIC1D_BEGIN_NAMESPACE
class BField;
class Current;

class EField : public Grid<CartVector> {
    Grid<CartVector>  E_prev;
    Grid<CovarVector> Bcovar;

public:
    ParamSet const params;
    Geometry const geomtr;

    explicit EField(ParamSet const &);
    EField &operator=(EField const &) = delete;

    [[nodiscard]] auto &grid_whole_domain_extent() const noexcept { return params.half_grid_whole_domain_extent; }
    [[nodiscard]] auto &grid_subdomain_extent() const noexcept { return params.half_grid_subdomain_extent; }

    void update(BField const &bfield, Current const &current, Real dt) noexcept;

private:
    void        mask(EField &, MaskingFunction const &masking_function) const;
    void        impl_update(EField &E_cart, Grid<CovarVector> const &Bcovar, Real cdtOsqrtg, Current const &J_cart, Real dt) const noexcept;
    static auto cart_to_covar(Grid<CovarVector> &Bcovar, BField const &B_cart) -> Grid<CovarVector> &;

    // attribute export facility
    //
    friend auto operator<<(hdf5::Group &obj, EField const &efield) -> decltype(obj);
    friend auto operator<<(hdf5::Dataset &obj, EField const &efield) -> decltype(obj);
    friend auto operator<<(hdf5::Group &&obj, EField const &efield) -> decltype(obj)
    {
        return std::move(obj << efield);
    }
    friend auto operator<<(hdf5::Dataset &&obj, EField const &efield) -> decltype(obj)
    {
        return std::move(obj << efield);
    }
};
PIC1D_END_NAMESPACE
