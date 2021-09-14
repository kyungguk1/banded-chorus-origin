/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "Recorder.h"

#include <limits>
#include <stdexcept>

HYBRID1D_BEGIN_NAMESPACE
namespace {
constexpr long large_int = std::numeric_limits<unsigned>::max();
}

Recorder::Recorder(unsigned const recording_frequency, parallel::mpi::Comm _comm)
: recording_frequency{ recording_frequency ? recording_frequency * Input::inner_Nt : large_int }
, comm{ std::move(_comm) }
{
    if (!comm->operator bool())
        throw std::invalid_argument{ std::string{ __PRETTY_FUNCTION__ } + " - invalid mpi::Comm" };
}

auto Recorder::get_space(std::vector<Scalar> const &payload) -> std::pair<hdf5::Space, hdf5::Space>
{
    static_assert(sizeof(Scalar) % sizeof(Real) == 0);
    static_assert(sizeof(Scalar) / sizeof(Real) == 1);

    auto mspace = hdf5::Space::simple(payload.size());
    mspace.select_all();

    auto fspace = hdf5::Space::simple(payload.size());
    fspace.select_all();

    return std::make_pair(mspace, fspace);
}
auto Recorder::get_space(std::vector<Vector> const &payload) -> std::pair<hdf5::Space, hdf5::Space>
{
    constexpr auto size = 3U;
    static_assert(sizeof(Vector) % sizeof(Real) == 0);
    static_assert(sizeof(Vector) / sizeof(Real) >= size);

    auto mspace = hdf5::Space::simple({ payload.size(), sizeof(Vector) / sizeof(Real) });
    mspace.select(H5S_SELECT_SET, { 0U, 0U }, { payload.size(), size });

    auto fspace = hdf5::Space::simple({ payload.size(), size });
    fspace.select_all();

    return std::make_pair(mspace, fspace);
}
auto Recorder::get_space(std::vector<Tensor> const &payload) -> std::pair<hdf5::Space, hdf5::Space>
{
    static_assert(sizeof(Tensor) % sizeof(Real) == 0);
    static_assert(sizeof(Tensor) / sizeof(Real) == 6);

    auto mspace = hdf5::Space::simple({ payload.size(), sizeof(Tensor) / sizeof(Real) });
    // diagonal
    mspace.select(H5S_SELECT_SET, { 0U, 0U }, { payload.size(), 3U });
    // off-diag
    // mspace.select(H5S_SELECT_OR, { 0U, 3U }, { payload.size(), 3U });

    auto fspace = hdf5::Space::simple({ payload.size(), 3U });
    fspace.select_all();

    return std::make_pair(mspace, fspace);
}
auto Recorder::get_space(std::vector<Particle::PSD> const &payload)
    -> std::pair<hdf5::Space, hdf5::Space>
{
    constexpr auto size = 3U;
    static_assert(sizeof(Particle::PSD) % sizeof(Real) == 0);
    static_assert(sizeof(Particle::PSD) / sizeof(Real) == size);

    auto mspace = hdf5::Space::simple({ payload.size(), size });
    mspace.select_all();

    auto fspace = hdf5::Space::simple({ payload.size(), size });
    fspace.select_all();

    return std::make_pair(mspace, fspace);
}
HYBRID1D_END_NAMESPACE
