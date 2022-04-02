/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC/VDF.h>
#include <algorithm>
#include <iterator>

LIBPIC_BEGIN_NAMESPACE
/// Test particle VDF
/// \details The sole job of this object is to dispense particles initialized with
///          the velocity and position passed by a TestParticleDesc object.
///
///          When all particles are exhausted, any further inquiry to the emit functions
///          returns a default-constructed Particle object which callers should ignore.
class TestParticleVDF : public VDF<TestParticleVDF> {
    friend VDF<TestParticleVDF>;

    static constexpr Real         m_Nrefcell_div_Ntotal{ 1 };
    KineticPlasmaDesc             desc;
    mutable std::vector<Particle> particles; // holder for remaining particles
    unsigned                      initial_number_of_test_particles;

public:
    /// Construct a test particle generator
    /// \tparam N The number of test particles.
    /// \param desc A TestParticleDesc object.
    /// \param geo A geometry object.
    /// \param domain_extent Spatial domain extent.
    /// \param c Light speed. (Not used here.)
    template <unsigned N>
    TestParticleVDF(TestParticleDesc<N> const &desc, Geometry const &geo, Range const &domain_extent, [[maybe_unused]] Real c)
    : VDF{ geo, domain_extent }, desc{ desc }, particles(N), initial_number_of_test_particles{ N }
    {
        std::transform(
            begin(desc.vel), end(desc.vel), begin(desc.pos), rbegin(particles) /*reverse order*/,
            [](Vector const &vel, CurviCoord const &pos) -> Particle {
                return { vel, pos };
            });
    }

private:
    [[nodiscard]] decltype(auto) impl_plasma_desc() const noexcept { return (this->desc); }

    [[nodiscard]] Scalar impl_n(CurviCoord const &) const { return 0; }
    [[nodiscard]] Vector impl_nV(CurviCoord const &) const { return { 0, 0, 0 }; }
    [[nodiscard]] Tensor impl_nvv(CurviCoord const &) const { return { 0, 0, 0, 0, 0, 0 }; }

    [[nodiscard]] Real impl_f0(Particle const &) const { return 0; }

    [[nodiscard]] std::vector<Particle> impl_emit(unsigned long) const;
    [[nodiscard]] Particle              impl_emit() const;
    [[nodiscard]] Particle              load() const;

public:
    // equilibrium physical distribution function
    //
    [[nodiscard]] Real f0(Vector const &, CurviCoord const &) const noexcept { return 0; }
    [[nodiscard]] Real f0(Particle const &ptl) const noexcept { return f0(ptl.vel, ptl.pos); }

    // marker particle distribution function
    //
    [[nodiscard]] Real g0(Vector const &, CurviCoord const &) const noexcept { return 1; }
    [[nodiscard]] Real g0(Particle const &ptl) const noexcept { return g0(ptl.vel, ptl.pos); }
};
LIBPIC_END_NAMESPACE
