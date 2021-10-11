/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include "Macros.h"
#include "Predefined.h"
#include <PIC/Grid.h>
#include <PIC/PlasmaDesc.h>
#include <PIC/Predefined.h>
#include <PIC/Range.h>
#include <PIC/Scalar.h>
#include <PIC/Tensor.h>
#include <PIC/Vector.h>

#include <array>
#include <cmath>
#include <tuple>
#include <type_traits>
#include <utility>

PIC1D_BEGIN_NAMESPACE
// input parameter import
//
#include <Inputs.h>

// grid definitions
//
using ScalarGrid = Grid<Scalar, Input::Nx / Input::number_of_subdomains, Pad>;
using VectorGrid = Grid<Vector, Input::Nx / Input::number_of_subdomains, Pad>;
using TensorGrid = Grid<Tensor, Input::Nx / Input::number_of_subdomains, Pad>;

// MARK:- Input Checks
//
namespace {
template <class Pred, class T, unsigned long... Is>
[[nodiscard]] constexpr auto is_all(
    Pred pred, std::array<T, sizeof...(Is)> A, std::index_sequence<Is...>) noexcept(noexcept(pred(std::declval<T>())))
    -> std::enable_if_t<std::is_invocable_r_v<bool, Pred, T const &>, bool>
{
    return (... && pred(std::get<Is>(A)));
}
template <class Pred, class T, unsigned long N>
[[nodiscard]] constexpr auto is_all(Pred pred, std::array<T, N> A) noexcept(noexcept(pred(std::declval<T>())))
    -> std::enable_if_t<std::is_invocable_r_v<bool, Pred, T const &>, bool>
{
    return is_all(pred, A, std::make_index_sequence<N>{});
}

template <class T, unsigned long N>
[[nodiscard]] constexpr bool is_all_positive(std::array<T, N> A)
{
    return is_all(
        [](T const &x) noexcept {
            return x > 0;
        },
        A);
}
template <class T, unsigned long N>
[[nodiscard]] constexpr bool is_all_nonnegative(std::array<T, N> A)
{
    return is_all(
        [](T const &x) noexcept {
            return x >= 0;
        },
        A);
}
template <class T, unsigned long N>
[[nodiscard]] constexpr bool is_all_non_zero(std::array<T, N> A)
{
    return is_all(
        [](T const &x) noexcept {
            return x != 0;
        },
        A);
}

template <class... Ts, class Int, Int... Is>
[[nodiscard]] constexpr bool check_Nc(std::tuple<Ts...> const &descs,
                                      std::integer_sequence<Int, Is...>) noexcept
{
    return is_all(
        [Nx = Input::Nx, denom = Input::number_of_worker_threads + 1](long const &x) noexcept {
            return x * Nx % denom == 0;
        },
        std::array<long, sizeof...(Ts)>{ std::get<Is>(descs).Nc... });
}
template <class... Ts>
[[nodiscard]] constexpr bool check_Nc(std::tuple<Ts...> const &descs) noexcept
{
    return check_Nc(descs, std::index_sequence_for<Ts...>{});
}
template <class... Ts, class Int, Int... Is>
[[nodiscard]] constexpr bool check_shape(std::tuple<Ts...> const &descs,
                                         std::integer_sequence<Int, Is...>) noexcept
{
    return is_all(
        [pad = Pad](ShapeOrder const &order) noexcept {
            return pad >= order;
        },
        std::array<ShapeOrder, sizeof...(Ts)>{ std::get<Is>(descs).shape_order... });
}
template <class... Ts>
[[nodiscard]] constexpr bool check_shape(std::tuple<Ts...> const &descs) noexcept
{
    return check_shape(descs, std::index_sequence_for<Ts...>{});
}
} // namespace

static_assert(Input::number_of_subdomains > 0, "number_of_subdomains should be a positive number");
static_assert(Input::number_of_distributed_particle_subdomain_clones > 0, "number_of_distributed_particle_subdomain_clones should be a positive number");
static_assert((1 + Input::number_of_worker_threads) % (Input::number_of_subdomains * Input::number_of_distributed_particle_subdomain_clones) == 0,
              "(1 + number_of_worker_threads) should be divisible by number_of_subdomains");

static_assert(Input::c > 0, "speed of light should be a positive number");
static_assert(Input::O0 > 0, "equatorial background magnetic field should be a positive number");
static_assert(Input::xi >= 0, "field inhomogeneity parameter should be a non-negative number");
static_assert(Input::Dx > 0, "grid size should be a positive number");
static_assert(Input::Nx > 0, "there should be at least 1 grid point");
static_assert(Input::Nx % Input::number_of_subdomains == 0, "Nx should be divisible by number_of_subdomains");
static_assert(Input::dt > 0, "time step should be a positive number");
static_assert(Input::inner_Nt > 0, "inner loop count should be a positive number");

static_assert(check_shape(Input::part_descs), "shape order should be less than or equal to the number of ghost cells");
PIC1D_END_NAMESPACE
