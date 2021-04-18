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

#include "SubdomainDelegate.h"

#include <random>
#include <stdexcept>
#include <utility>

H1D::SubdomainDelegate::message_dispatch_t H1D::SubdomainDelegate::dispatch{
    H1D::ParamSet::number_of_subdomains};
H1D::SubdomainDelegate::SubdomainDelegate(unsigned const rank, unsigned const size)
: comm{dispatch.comm(rank)}
, size{size}
, left_{(size + rank - 1) % size}
, right{(size + rank + 1) % size}
{
    if (size > ParamSet::number_of_subdomains)
        throw std::invalid_argument{__PRETTY_FUNCTION__};
}

// MARK: Interface
//
void H1D::SubdomainDelegate::once(Domain &domain) const
{
    std::mt19937                     g{123 + static_cast<unsigned>(comm.rank)};
    std::uniform_real_distribution<> d{-1, 1};
    for (Vector &v : domain.bfield) {
        v.x += d(g) * Debug::initial_bfield_noise_amplitude;
        v.y += d(g) * Debug::initial_bfield_noise_amplitude;
        v.z += d(g) * Debug::initial_bfield_noise_amplitude;
    }
}

void H1D::SubdomainDelegate::pass(Domain const &domain, PartBucket &L_bucket,
                                  PartBucket &R_bucket) const
{
    // pass across boundaries
    //
    {
        auto tk1 = comm.send(std::move(L_bucket), left_);
        auto tk2 = comm.send(std::move(R_bucket), right);
        L_bucket = comm.recv<PartBucket>(right);
        R_bucket = comm.recv<PartBucket>(left_);
        // std::move(tk1).wait();
        // std::move(tk2).wait();
    }

    // adjust coordinates
    //
    Delegate::pass(domain, L_bucket, R_bucket);
}
void H1D::SubdomainDelegate::pass(Domain const &, ColdSpecies &sp) const
{
    pass(sp.mom0_full);
    pass(sp.mom1_full);
}
void H1D::SubdomainDelegate::pass(Domain const &, BField &bfield) const
{
    if constexpr (Debug::zero_out_electromagnetic_field) {
        bfield.fill(bfield.geomtr.B0);
    }
    pass(bfield);
}
void H1D::SubdomainDelegate::pass(Domain const &, EField &efield) const
{
    if constexpr (Debug::zero_out_electromagnetic_field) {
        efield.fill(Vector{});
    }
    pass(efield);
}
void H1D::SubdomainDelegate::pass(Domain const &, Charge &charge) const
{
    pass(charge);
}
void H1D::SubdomainDelegate::pass(Domain const &, Current &current) const
{
    pass(current);
}
void H1D::SubdomainDelegate::gather(Domain const &, Charge &charge) const
{
    gather(charge);
}
void H1D::SubdomainDelegate::gather(Domain const &, Current &current) const
{
    gather(current);
}
void H1D::SubdomainDelegate::gather(Domain const &, PartSpecies &sp) const
{
    gather(sp.moment<0>());
    gather(sp.moment<1>());
    gather(sp.moment<2>());
}

template <class T, long N> void H1D::SubdomainDelegate::pass(GridQ<T, N> &grid) const
{
    // from inside out
    //
    auto tk1 = comm.send<T const *>(grid.begin(), left_);
    auto tk2 = comm.send<T const *>(grid.end(), right);
    //
    comm.recv<T const *>(right).unpack(
        [](T const *right, T *last) {
            for (long i = 0; i < Pad; ++i) {
                last[i] = right[i];
            }
        },
        grid.end());
    //
    comm.recv<T const *>(left_).unpack(
        [](T const *left, T *first) {
            for (long i = -1; i >= -Pad; --i) {
                first[i] = left[i];
            }
        },
        grid.begin());
    //
    std::move(tk1).wait();
    std::move(tk2).wait();
}
template <class T, long N> void H1D::SubdomainDelegate::gather(GridQ<T, N> &grid) const
{
    // from outside in
    //
    auto tk1 = comm.send<T const *>(grid.begin(), left_);
    auto tk2 = comm.send<T const *>(grid.end(), right);
    //
    comm.recv<T const *>(right).unpack(
        [](T const *right, T *last) {
            for (long i = -Pad; i < 0; ++i) {
                last[i] += right[i];
            }
        },
        grid.end());
    //
    comm.recv<T const *>(left_).unpack(
        [](T const *left, T *first) {
            for (long i = Pad - 1; i >= 0; --i) {
                first[i] += left[i];
            }
        },
        grid.begin());
    //
    std::move(tk1).wait();
    std::move(tk2).wait();
}
