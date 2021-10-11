/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "TestParticleVDF.h"

LIBPIC_BEGIN_NAMESPACE
auto TestParticleVDF::impl_emit(unsigned long const n) const -> std::vector<Particle>
{
    std::vector<Particle> ptls(n);
    for (auto &ptl : ptls)
        ptl = emit();
    return ptls;
}
auto TestParticleVDF::impl_emit() const -> Particle
{
    Particle ptl = load();
    {
        ptl.psd = { 1, -1, -1 };
    }
    return ptl;
}
auto TestParticleVDF::load() const -> Particle
{
    if (particles.empty())
        return {}; // this assumes that the default-constructed particle object is not consumed by callers

    auto ptl = particles.back();
    particles.pop_back();
    return ptl;
}
LIBPIC_END_NAMESPACE
