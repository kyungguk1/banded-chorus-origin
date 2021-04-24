/*
 * Copyright (c) 2019-2021, Kyungguk Min
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

// MARK:- P1D::SubdomainDelegate
//
P1D::SubdomainDelegate::SubdomainDelegate(parallel::mpi::Comm _comm) : comm{std::move(_comm), tag}
{
    if (!comm->operator bool())
        throw std::invalid_argument{__PRETTY_FUNCTION__};

    int const size = comm.size();
    int const rank = comm->rank();
    left_          = rank_t{(size + rank - 1) % size};
    right          = rank_t{(size + rank + 1) % size};
}

void P1D::SubdomainDelegate::once(Domain &domain) const
{
    std::mt19937                     g{494983U + static_cast<unsigned>(comm->rank())};
    std::uniform_real_distribution<> d{-1, 1};
    for (Vector &v : domain.efield) {
        v.x += d(g) * Debug::initial_efield_noise_amplitude;
        v.y += d(g) * Debug::initial_efield_noise_amplitude;
        v.z += d(g) * Debug::initial_efield_noise_amplitude;
    }
}

void P1D::SubdomainDelegate::pass(Domain const &, ColdSpecies &sp) const
{
    pass(sp.mom0_full);
    pass(sp.mom1_full);
}
void P1D::SubdomainDelegate::pass(Domain const &, BField &bfield) const
{
    if constexpr (Debug::zero_out_electromagnetic_field) {
        bfield.fill(bfield.geomtr.B0);
    } else if constexpr (Input::is_electrostatic) { // zero-out transverse components
        for (Vector &v : bfield) {
            v.y = bfield.geomtr.B0.y;
            v.z = bfield.geomtr.B0.z;
        }
    }
    pass(bfield);
}
void P1D::SubdomainDelegate::pass(Domain const &, EField &efield) const
{
    if constexpr (Debug::zero_out_electromagnetic_field) {
        efield.fill(Vector{});
    } else if constexpr (Input::is_electrostatic) { // zero-out transverse components
        for (Vector &v : efield) {
            v.y = v.z = 0;
        }
    }
    pass(efield);
}
void P1D::SubdomainDelegate::pass(Domain const &, Current &current) const
{
    pass(current);
}
void P1D::SubdomainDelegate::gather(Domain const &, Current &current) const
{
    gather(current);
}
void P1D::SubdomainDelegate::gather(Domain const &, PartSpecies &sp) const
{
    gather(sp.moment<0>());
    gather(sp.moment<1>());
    gather(sp.moment<2>());
}
void P1D::SubdomainDelegate::pass(Domain const &domain, PartBucket &L_bucket, PartBucket &R_bucket) const
{
    // pass across boundaries
    //
    {
        auto tk1 = comm.ibsend(std::move(L_bucket), left_);
        auto tk2 = comm.ibsend(std::move(R_bucket), right);
        L_bucket = comm.recv<3>({}, right);
        R_bucket = comm.recv<3>({}, left_);
        std::move(tk1).wait();
        std::move(tk2).wait();
    }

    // adjust coordinates
    //
    Delegate::pass(domain, L_bucket, R_bucket);
}
template <class T, long N> void P1D::SubdomainDelegate::pass(GridQ<T, N> &grid) const
{
    // from inside out
    //
    for (long b = 0, e = -1; b < Pad; ++b, --e) {
        {
            auto tk       = comm.issend<T>(grid.begin()[b], left_);
            grid.end()[b] = comm.recv<T>(right);
            std::move(tk).wait();
        }
        {
            auto tk         = comm.issend<T>(grid.end()[e], right);
            grid.begin()[e] = comm.recv<T>(left_);
            std::move(tk).wait();
        }
    }
}
template <class T, long N> void P1D::SubdomainDelegate::gather(GridQ<T, N> &grid) const
{
    // from outside in
    //
    for (long b = -Pad, e = Pad - 1; b < 0; ++b, --e) {
        {
            auto tk = comm.issend<T>(grid.begin()[b], left_);
            grid.end()[b] += comm.recv<T>(right);
            std::move(tk).wait();
        }
        {
            auto tk = comm.issend<T>(grid.end()[e], right);
            grid.begin()[e] += comm.recv<T>(left_);
            std::move(tk).wait();
        }
    }
}
