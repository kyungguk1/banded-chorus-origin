/*
 * Copyright (c) 2021-2022, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "catch2/catch.hpp"

#define LIBPIC_INLINE_VERSION 1
//#include <PIC/RelativisticVDF/LossconeVDF.h>
#include <PIC/RelativisticVDF/MaxwellianVDF.h>
//#include <PIC/RelativisticVDF/PartialShellVDF.h>
#include <PIC/RelativisticVDF/TestParticleVDF.h>
#include <PIC/UTL/println.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <numeric>

namespace {
constexpr bool dump_samples         = false;
constexpr bool enable_moment_checks = false;

using Particle = RelativisticParticle;

[[nodiscard]] bool operator==(Scalar const &a, Scalar const &b) noexcept
{
    return *a == Approx{ *b }.margin(1e-15);
}
template <class T, class U>
[[nodiscard]] bool operator==(Detail::VectorTemplate<T, double> const &a, Detail::VectorTemplate<U, double> const &b) noexcept
{
    return a.x == Approx{ b.x }.margin(1e-15)
        && a.y == Approx{ b.y }.margin(1e-15)
        && a.z == Approx{ b.z }.margin(1e-15);
}
template <class T1, class T2, class U1, class U2>
[[nodiscard]] bool operator==(Detail::TensorTemplate<T1, T2> const &a, Detail::TensorTemplate<U1, U2> const &b) noexcept
{
    return a.lo() == b.lo() && a.hi() == b.hi();
}
template <class T1, class T2, class U1, class U2>
[[nodiscard]] bool operator==(Detail::FourVectorTemplate<T1, T2> const &a, Detail::FourVectorTemplate<U1, U2> const &b) noexcept
{
    return a.t == Approx{ b.t }.margin(1e-15) && a.s == b.s;
}
template <class T1, class T2, class T3, class U1, class U2, class U3>
[[nodiscard]] bool operator==(Detail::FourTensorTemplate<T1, T2, T3> const &a, Detail::FourTensorTemplate<U1, U2, U3> const &b) noexcept
{
    return a.tt == Approx{ b.tt }.margin(1e-15)
        && a.ts == b.ts
        && a.ss == b.ss;
}
[[nodiscard]] bool operator==(CurviCoord const &a, CurviCoord const &b) noexcept
{
    return a.q1 == b.q1;
}
} // namespace
using ::operator==;

TEST_CASE("Test LibPIC::RelativisticVDF::TestParticleVDF", "[LibPIC::RelativisticVDF::TestParticleVDF]")
{
    Real const O0 = 1., op = 10 * O0, c = op;
    Real const xi = 0, D1 = 1;
    long const q1min = -7, q1max = 15;
    auto const geo   = Geometry{ xi, D1, O0 };
    auto const Nptls = 2;
    auto const desc  = TestParticleDesc<Nptls>{
        { -O0, op },
        { MFAVector{ 1, 2, 3 }, { 3, 4, 5 } },
        { CurviCoord{ q1min }, CurviCoord{ q1max } }
    };
    auto const vdf = RelativisticTestParticleVDF(desc, geo, { q1min, q1max - q1min }, c);

    CHECK(serialize(desc) == serialize(vdf.plasma_desc()));

    CHECK(vdf.initial_number_of_test_particles == Nptls);
    CHECK(vdf.Nrefcell_div_Ntotal() == Approx{ 1.0 }.epsilon(1e-10));

    for (long q1 = q1min; q1 <= q1max; ++q1) {
        CurviCoord const pos(q1);

        CHECK(vdf.n0(pos) == Scalar{ 0 });
        CHECK(geo.cart_to_mfa(vdf.nV0(pos), pos) == MFAVector{ 0, 0, 0 });
        CHECK(geo.cart_to_mfa(vdf.nuv0(pos), pos) == FourMFATensor{ 0, { 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 } });
    }

    // sampling
    auto const n_samples = Nptls + 1;
    auto const particles = vdf.emit(n_samples);

    for (unsigned i = 0; i < Nptls; ++i) {
        auto const &ptl = particles[i];

        REQUIRE(ptl.psd.weight == 0);
        REQUIRE(ptl.psd.real_f == 0);
        REQUIRE(ptl.psd.marker == 1);

        REQUIRE(ptl.pos == desc.pos[i]);
        REQUIRE(geo.cart_to_mfa(ptl.gcgvel, ptl.pos) == lorentz_boost<-1>(FourMFAVector{ c, {} }, desc.vel[i] / c, Real{ ptl.gcgvel.t } / c));

        REQUIRE(vdf.real_f0(ptl) == 0);
        REQUIRE(vdf.g0(ptl) == 1);
    }
    {
        auto const &ptl = particles.back();

        REQUIRE(std::isnan(*ptl.gcgvel.t));
        REQUIRE(std::isnan(ptl.gcgvel.s.x));
        REQUIRE(std::isnan(ptl.gcgvel.s.y));
        REQUIRE(std::isnan(ptl.gcgvel.s.z));
        REQUIRE(std::isnan(ptl.pos.q1));
    }

    if constexpr (dump_samples) {
        static_assert(n_samples > 0);
        std::ofstream os{ "/Users/kyungguk/Downloads/RelativisticTestParticleVDF.m" };
        os.setf(os.fixed);
        os.precision(20);
        println(os, '{');
        std::for_each(begin(particles), std::prev(end(particles)), [&os](auto const &ptl) {
            println(os, "    ", ptl, ", ");
        });
        println(os, "    ", particles.back());
        println(os, '}');
        os.close();
    }
}

TEST_CASE("Test LibPIC::RelativisticVDF::MaxwellianVDF::Homogeneous", "[LibPIC::RelativisticVDF::MaxwellianVDF::Homogeneous]")
{
    Real const O0 = 1., op = 4 * O0, c = op, beta1 = 1.5, T2OT1 = 5.35, Vd = 0;
    Real const xi = 0, D1 = 1.87, psd_refresh_frequency = 0;
    long const q1min = -7, q1max = 15;
    auto const geo     = Geometry{ xi, D1, O0 };
    auto const kinetic = KineticPlasmaDesc({ -O0, op }, 10, ShapeOrder::CIC, psd_refresh_frequency, ParticleScheme::delta_f, 0, 2.1);
    auto const desc    = BiMaxPlasmaDesc(kinetic, beta1, T2OT1);
    auto const vdf     = RelativisticMaxwellianVDF(desc, geo, { q1min, q1max - q1min }, c);

    auto const g_desc = BiMaxPlasmaDesc({ { -O0, op }, 10, ShapeOrder::CIC }, beta1 * desc.marker_temp_ratio, T2OT1);
    auto const g_vdf  = RelativisticMaxwellianVDF(g_desc, geo, { q1min, q1max - q1min }, c);

    CHECK(serialize(desc) == serialize(vdf.plasma_desc()));

    // check equilibrium macro variables
    CHECK(vdf.Nrefcell_div_Ntotal() == Approx{ 1.0 / (q1max - q1min) }.epsilon(1e-10));

    for (long q1 = q1min; q1 <= q1max; ++q1) {
        CurviCoord const pos(q1);

        auto const n0_ref   = 1;
        auto const nV0_ref  = n0_ref * Vector{ Vd, 0, 0 };
        auto const nuv0_ref = FourTensor{
            19.68844169847468,
            { 0, 0, 0 },
            { 0.6024514298876814, 2.913078061288704, 2.913078061288704, 0, 0, 0 },
        };

        REQUIRE(vdf.n0(pos) == Scalar{ n0_ref });
        REQUIRE(geo.cart_to_mfa(vdf.nV0(pos), pos) == nV0_ref);
        REQUIRE(geo.cart_to_mfa(vdf.nuv0(pos), pos) == nuv0_ref);
    }

    // sampling
    auto const n_samples = 100000U;
    auto       particles = vdf.emit(n_samples);

    // moments
    if constexpr (enable_moment_checks) {
        auto const particle_density = std::accumulate(
            begin(particles), end(particles), Real{}, [&](Real const sum, Particle const &ptl) {
                return sum + 1 * (ptl.psd.real_f / ptl.psd.marker + desc.scheme * ptl.psd.weight) / n_samples;
            });
        CHECK(particle_density == Approx{ 1 }.epsilon(1e-2));

        auto const particle_flux = std::accumulate(
            begin(particles), end(particles), CartVector{}, [&](CartVector const &sum, Particle const &ptl) {
                auto const sqrt_eta = std::sqrt(geo.Bmag_div_B0(ptl.pos) / (1 - desc.T2_T1 * (1 - geo.Bmag_div_B0(ptl.pos))));
                auto const vel      = ptl.velocity(c) / CartVector{ 1, sqrt_eta, sqrt_eta } / std::sqrt(desc.beta1);
                return sum + vel * (ptl.psd.real_f / ptl.psd.marker + desc.scheme * ptl.psd.weight) / n_samples;
            });
        CHECK(particle_flux.x == Approx{ 0 }.margin(1e-2));
        CHECK(particle_flux.y == Approx{ 0 }.margin(1e-2));
        CHECK(particle_flux.z == Approx{ 0 }.margin(1e-2));

        auto const stress_energy = std::accumulate(
            begin(particles), end(particles), FourCartTensor{}, [&](FourCartTensor const &sum, Particle const &ptl) {
                auto const sqrt_eta = std::sqrt(geo.Bmag_div_B0(ptl.pos) / (1 - desc.T2_T1 * (1 - geo.Bmag_div_B0(ptl.pos))));
                auto const gcgv     = ptl.gcgvel / FourCartVector{ 1, CartVector{ 1, sqrt_eta, sqrt_eta } };
                auto const v        = gcgv.s / Real{ gcgv.t } * c;
                auto const u        = gcgv.s;
                auto const ft       = FourCartTensor{
                    c * gcgv.t / (c * c),
                    c * gcgv.s / c / std::sqrt(desc.beta1),
                    CartTensor{ v.x * u.x, v.y * u.y, v.z * u.z, v.x * u.y, v.y * u.z, v.z * u.x } / desc.beta1
                };
                return sum + ft * (ptl.psd.real_f / ptl.psd.marker + desc.scheme * ptl.psd.weight) / n_samples;
            });
        auto const nuv0 = vdf.nuv0({}) / FourCartTensor{ c * c, CartVector{ std::sqrt(desc.beta1) }, CartTensor{ desc.beta1 } };
        CHECK(*stress_energy.tt == Approx{ *nuv0.tt }.epsilon(1e-2));
        CHECK(stress_energy.ts.x == Approx{ nuv0.ts.x }.margin(1e-2));
        CHECK(stress_energy.ts.y == Approx{ nuv0.ts.y }.margin(1e-2));
        CHECK(stress_energy.ts.z == Approx{ nuv0.ts.z }.margin(1e-2));
        CHECK(stress_energy.ss.xx == Approx{ nuv0.ss.xx }.epsilon(1e-2));
        CHECK(stress_energy.ss.yy == Approx{ nuv0.ss.yy }.epsilon(1e-2));
        CHECK(stress_energy.ss.zz == Approx{ nuv0.ss.zz }.epsilon(1e-2));
        CHECK(stress_energy.ss.xy == Approx{ nuv0.ss.xy }.margin(1e-2));
        CHECK(stress_energy.ss.yz == Approx{ nuv0.ss.yz }.margin(1e-2));
        CHECK(stress_energy.ss.zx == Approx{ nuv0.ss.zx }.margin(1e-2));
    }

    static_assert(n_samples > 100);
    for (unsigned long i = 0; i < 100; ++i) {
        Particle const &ptl = particles[i];

        REQUIRE(ptl.psd.weight == Approx{ desc.initial_weight * desc.scheme + (1 - desc.scheme) * vdf.f0(ptl) / g_vdf.f0(ptl) }.margin(1e-15));
        REQUIRE(ptl.psd.marker == Approx{ g_vdf.f0(ptl) }.epsilon(1e-14));
        REQUIRE(ptl.psd.real_f == Approx{ vdf.f0(ptl) * desc.scheme + ptl.psd.weight * ptl.psd.marker }.epsilon(1e-10));
        REQUIRE(vdf.real_f0(ptl) == Approx{ vdf.f0(ptl) }.epsilon(1e-14));
        REQUIRE(ptl.gcgvel.t == Approx{ std::sqrt(c * c + dot(ptl.gcgvel.s, ptl.gcgvel.s)) }.epsilon(1e-10));

        auto const gd   = c / std::sqrt((c - Vd) * (c + Vd));
        auto const n0   = *vdf.n0(ptl.pos) / gd;
        auto const vth1 = std::sqrt(beta1);
        auto const vth2 = vth1 * std::sqrt(T2OT1);
        auto const u_co = lorentz_boost<+1>(geo.cart_to_mfa(ptl.gcgvel, ptl.pos), Vd / c, gd).s;
        auto const u1   = u_co.x;
        auto const u2   = std::sqrt(u_co.y * u_co.y + u_co.z * u_co.z);
        auto const f_ref
            = n0 * std::exp(-u1 * u1 / (vth1 * vth1) - u2 * u2 / (vth2 * vth2))
            / (4 * M_PI_2 / M_2_SQRTPI * vth1 * vth2 * vth2);
        auto const g_ref
            = n0 * std::exp(-u1 * u1 / (vth1 * vth1 * desc.marker_temp_ratio) - u2 * u2 / (vth2 * vth2 * desc.marker_temp_ratio))
            / (4 * M_PI_2 / M_2_SQRTPI * vth1 * vth2 * vth2 * desc.marker_temp_ratio * std::sqrt(desc.marker_temp_ratio));

        REQUIRE(vdf.f0(ptl) == Approx{ f_ref }.epsilon(1e-10));
        REQUIRE(vdf.g0(ptl) == Approx{ g_ref }.epsilon(1e-10));
    }

    if constexpr (dump_samples) {
        static_assert(n_samples > 0);
        std::ofstream os{ "/Users/kyungguk/Downloads/RelativisticMaxwellianVDF-homogeneous.m" };
        os.setf(os.fixed);
        os.precision(20);
        println(os, '{');
        for (unsigned long i = 0; i < particles.size() - 1; ++i) {
            println(os, "    ", particles[i], ", ");
        }
        println(os, "    ", particles.back());
        println(os, '}');
        os.close();
    }
}

TEST_CASE("Test LibPIC::RelativisticVDF::MaxwellianVDF::Inhomogeneous", "[LibPIC::RelativisticVDF::MaxwellianVDF::Inhomogeneous]")
{
    Real const O0 = 1., op = 4 * O0, c = op, beta1_eq = 1.5, T2OT1_eq = 5.35;
    Real const xi = .876, xiD1q1max = M_PI_2 * 0.8, psd_refresh_frequency = 0;
    long const q1min = -7, q1max = 15;
    auto const D1      = xiD1q1max / (xi * std::max(std::abs(q1min), std::abs(q1max)));
    auto const geo     = Geometry{ xi, D1, O0 };
    auto const kinetic = KineticPlasmaDesc({ -O0, op }, 10, ShapeOrder::CIC, psd_refresh_frequency, ParticleScheme::delta_f, 0, 2.1);
    auto const desc    = BiMaxPlasmaDesc(kinetic, beta1_eq, T2OT1_eq);
    auto const vdf     = RelativisticMaxwellianVDF(desc, geo, { q1min, q1max - q1min }, c);

    auto const g_desc = BiMaxPlasmaDesc({ { -O0, op }, 10, ShapeOrder::CIC }, beta1_eq * desc.marker_temp_ratio, T2OT1_eq);
    auto const g_vdf  = RelativisticMaxwellianVDF(g_desc, geo, { q1min, q1max - q1min }, c);

    CHECK(serialize(desc) == serialize(vdf.plasma_desc()));

    // check equilibrium macro variables
    CHECK(vdf.Nrefcell_div_Ntotal() == Approx{ 0.0796995979174554 }.epsilon(1e-10));

    std::array const etas{
        0.4287882239954878, 0.49761722933466956, 0.5815167468245395, 0.6800561885442317, 0.788001590917828,
        0.8920754659248894, 0.9704416860618577, 1., 0.9704416860618577, 0.8920754659248894, 0.788001590917828,
        0.6800561885442317, 0.5815167468245395, 0.49761722933466956, 0.4287882239954878, 0.3733660301900688,
        0.32911727304544663, 0.29391455144617107, 0.26596447719479277, 0.24383880547267717, 0.22643378958243887,
        0.21291237930064547, 0.2026501794964986
    };
    std::array<FourMFATensor, std::size(etas)> const nuv0s{
        FourMFATensor{ 7.6796439042989863566, { 0, 0, 0 }, { 0.27911539238104404737, 0.61594144886886081913, 0.61594144886886081913, 0, 0, 0 } },
        FourMFATensor{ 9.0259795339782940005, { 0, 0, 0 }, { 0.32050341039122298703, 0.81325808234008256647, 0.81325808234008256647, 0, 0, 0 } },
        FourMFATensor{ 10.706184755021030952, { 0, 0, 0 }, { 0.36993938519636432316, 1.0855481415092880226, 1.0855481415092880226, 0, 0, 0 } },
        FourMFATensor{ 12.732781177900619696, { 0, 0, 0 }, { 0.42668819651599032561, 1.4477412024058688989, 1.4477412024058688989, 0, 0, 0 } },
        FourMFATensor{ 15.016692616105265401, { 0, 0, 0 }, { 0.48735966029137123279, 1.8943437716588384934, 1.8943437716588384934, 0, 0, 0 } },
        FourMFATensor{ 17.279918470199586267, { 0, 0, 0 }, { 0.54449711801919620235, 2.3717632904336083399, 2.3717632904336083399, 0, 0, 0 } },
        FourMFATensor{ 19.022668266008288640, { 0, 0, 0 }, { 0.58670448769353888974, 2.7602655816622587714, 2.7602655816622587714, 0, 0, 0 } },
        FourMFATensor{ 19.688441698474679953, { 0, 0, 0 }, { 0.60245142988768141112, 2.9130780612887035019, 2.9130780612887035019, 0, 0, 0 } },
        FourMFATensor{ 19.022668266008288640, { 0, 0, 0 }, { 0.58670448769353888974, 2.7602655816622587714, 2.7602655816622587714, 0, 0, 0 } },
        FourMFATensor{ 17.279918470199586267, { 0, 0, 0 }, { 0.54449711801919620235, 2.3717632904336083399, 2.3717632904336083399, 0, 0, 0 } },
        FourMFATensor{ 15.016692616105265401, { 0, 0, 0 }, { 0.48735966029137123279, 1.8943437716588384934, 1.8943437716588384934, 0, 0, 0 } },
        FourMFATensor{ 12.732781177900619696, { 0, 0, 0 }, { 0.42668819651599032561, 1.4477412024058688989, 1.4477412024058688989, 0, 0, 0 } },
        FourMFATensor{ 10.706184755021030952, { 0, 0, 0 }, { 0.36993938519636432316, 1.0855481415092880226, 1.0855481415092880226, 0, 0, 0 } },
        FourMFATensor{ 9.0259795339782940005, { 0, 0, 0 }, { 0.32050341039122298703, 0.81325808234008256647, 0.81325808234008256647, 0, 0, 0 } },
        FourMFATensor{ 7.6796439042989863566, { 0, 0, 0 }, { 0.27911539238104404737, 0.61594144886886081913, 0.61594144886886081913, 0, 0, 0 } },
        FourMFATensor{ 6.6171230649494150455, { 0, 0, 0 }, { 0.24520396653403486731, 0.47492035609505728333, 0.47492035609505728333, 0, 0, 0 } },
        FourMFATensor{ 5.7829467934930152140, { 0, 0, 0 }, { 0.21773090191556496165, 0.37422770055057413829, 0.37422770055057413829, 0, 0, 0 } },
        FourMFATensor{ 5.1284442913639001205, { 0, 0, 0 }, { 0.19560829100507137746, 0.30192268264068905514, 0.30192268264068905514, 0, 0, 0 } },
        FourMFATensor{ 4.6146462260100493680, { 0, 0, 0 }, { 0.17786822147555930718, 0.24957915752025311429, 0.24957915752025311429, 0, 0, 0 } },
        FourMFATensor{ 4.2116440004148207876, { 0, 0, 0 }, { 0.16371093623219648561, 0.21139778865898520288, 0.21139778865898520288, 0, 0, 0 } },
        FourMFATensor{ 3.8969648050428511432, { 0, 0, 0 }, { 0.15250132374346736519, 0.18342295339037723023, 0.18342295339037723023, 0, 0, 0 } },
        FourMFATensor{ 3.6539349702175152323, { 0, 0, 0 }, { 0.14374758536644874352, 0.16296207621484651296, 0.16296207621484651296, 0, 0, 0 } },
        FourMFATensor{ 3.4703283142123186877, { 0, 0, 0 }, { 0.13707689292660604763, 0.14818475209338716203, 0.14818475209338716203, 0, 0, 0 } },
    };
    for (unsigned i = 0; i < std::size(etas); ++i) {
        auto const &eta = etas.at(i);
        auto const  q1  = q1min + i;
        auto const  pos = CurviCoord(q1);

        auto const n0_ref   = eta;
        auto const nV0_ref  = Vector{};
        auto const nuv0_ref = nuv0s.at(i);

        REQUIRE(vdf.n0(pos) == Scalar{ n0_ref });
        REQUIRE(geo.cart_to_mfa(vdf.nV0(pos), pos) == nV0_ref);
        REQUIRE(geo.cart_to_mfa(vdf.nuv0(pos), pos) == nuv0_ref);
    }

    // sampling
    auto const n_samples = 100000U;
    auto       particles = vdf.emit(n_samples);

    // moments
    if constexpr (enable_moment_checks) {
        auto const particle_density = std::accumulate(
            begin(particles), end(particles), Real{}, [&](Real const sum, Particle const &ptl) {
                return sum + 1 * (ptl.psd.real_f / ptl.psd.marker + desc.scheme * ptl.psd.weight) / n_samples;
            });
        CHECK(particle_density == Approx{ 1 }.epsilon(1e-2));

        auto const particle_flux = std::accumulate(
            begin(particles), end(particles), CartVector{}, [&](CartVector const &sum, Particle const &ptl) {
                auto const sqrt_eta = std::sqrt(geo.Bmag_div_B0(ptl.pos) / (1 - desc.T2_T1 * (1 - geo.Bmag_div_B0(ptl.pos))));
                auto const vel      = ptl.velocity(c) / CartVector{ 1, sqrt_eta, sqrt_eta } / std::sqrt(desc.beta1);
                return sum + vel * (ptl.psd.real_f / ptl.psd.marker + desc.scheme * ptl.psd.weight) / n_samples;
            });
        CHECK(particle_flux.x == Approx{ 0 }.margin(1e-2));
        CHECK(particle_flux.y == Approx{ 0 }.margin(1e-2));
        CHECK(particle_flux.z == Approx{ 0 }.margin(1e-2));

        auto const stress_energy = std::accumulate(
            begin(particles), end(particles), FourCartTensor{}, [&](FourCartTensor const &sum, Particle const &ptl) {
                auto const sqrt_eta = std::sqrt(geo.Bmag_div_B0(ptl.pos) / (1 - desc.T2_T1 * (1 - geo.Bmag_div_B0(ptl.pos))));
                auto const gcgv     = ptl.gcgvel / FourCartVector{ 1, CartVector{ 1, sqrt_eta, sqrt_eta } };
                auto const v        = gcgv.s / Real{ gcgv.t } * c;
                auto const u        = gcgv.s;
                auto const ft       = FourCartTensor{
                    c * gcgv.t / (c * c),
                    c * gcgv.s / c / std::sqrt(desc.beta1),
                    CartTensor{ v.x * u.x, v.y * u.y, v.z * u.z, v.x * u.y, v.y * u.z, v.z * u.x } / desc.beta1
                };
                return sum + ft * (ptl.psd.real_f / ptl.psd.marker + desc.scheme * ptl.psd.weight) / n_samples;
            });
        auto const nuv0 = vdf.nuv0({}) / FourCartTensor{ c * c, CartVector{ std::sqrt(desc.beta1) }, CartTensor{ desc.beta1 } };
        CHECK(*stress_energy.tt == Approx{ *nuv0.tt }.epsilon(1e-1));
        CHECK(stress_energy.ts.x == Approx{ nuv0.ts.x }.margin(1e-2));
        CHECK(stress_energy.ts.y == Approx{ nuv0.ts.y }.margin(1e-2));
        CHECK(stress_energy.ts.z == Approx{ nuv0.ts.z }.margin(1e-2));
        CHECK(stress_energy.ss.xx == Approx{ nuv0.ss.xx }.epsilon(1e-1));
        CHECK(stress_energy.ss.yy == Approx{ nuv0.ss.yy }.epsilon(1e-1));
        CHECK(stress_energy.ss.zz == Approx{ nuv0.ss.zz }.epsilon(1e-1));
        CHECK(stress_energy.ss.xy == Approx{ nuv0.ss.xy }.margin(1e-2));
        CHECK(stress_energy.ss.yz == Approx{ nuv0.ss.yz }.margin(1e-2));
        CHECK(stress_energy.ss.zx == Approx{ nuv0.ss.zx }.margin(1e-2));
    }

    static_assert(n_samples > 100);
    for (unsigned long i = 0; i < 100; ++i) {
        Particle const &ptl = particles[i];

        REQUIRE(ptl.psd.weight == Approx{ desc.initial_weight * desc.scheme + (1 - desc.scheme) * vdf.f0(ptl) / g_vdf.f0(ptl) }.margin(1e-15));
        REQUIRE(ptl.psd.marker == Approx{ g_vdf.f0(ptl) }.epsilon(1e-14));
        REQUIRE(ptl.psd.real_f == Approx{ vdf.f0(ptl) * desc.scheme + ptl.psd.weight * ptl.psd.marker }.epsilon(1e-10));
        REQUIRE(vdf.real_f0(ptl) == Approx{ vdf.f0(ptl) }.epsilon(1e-14));
        REQUIRE(ptl.gcgvel.t == Approx{ std::sqrt(c * c + dot(ptl.gcgvel.s, ptl.gcgvel.s)) }.epsilon(1e-10));

        auto const n0   = *vdf.n0(ptl.pos);
        auto const vth1 = std::sqrt(beta1_eq);
        auto const vth2 = vth1 * std::sqrt(T2OT1_eq * n0);
        auto const uel  = geo.cart_to_mfa(ptl.gcgvel, ptl.pos).s;
        auto const u1   = uel.x;
        auto const u2   = std::sqrt(uel.y * uel.y + uel.z * uel.z);
        auto const f_ref
            = n0 * std::exp(-u1 * u1 / (vth1 * vth1) - u2 * u2 / (vth2 * vth2))
            / (4 * M_PI_2 / M_2_SQRTPI * vth1 * vth2 * vth2);
        auto const g_ref
            = n0 * std::exp(-u1 * u1 / (vth1 * vth1 * desc.marker_temp_ratio) - u2 * u2 / (vth2 * vth2 * desc.marker_temp_ratio))
            / (4 * M_PI_2 / M_2_SQRTPI * vth1 * vth2 * vth2 * desc.marker_temp_ratio * std::sqrt(desc.marker_temp_ratio));

        REQUIRE(vdf.f0(ptl) == Approx{ f_ref }.epsilon(1e-10));
        REQUIRE(vdf.g0(ptl) == Approx{ g_ref }.epsilon(1e-10));
    }

    if constexpr (dump_samples) {
        std::ofstream os{ "/Users/kyungguk/Downloads/RelativisticMaxwellianVDF-inhomogeneous.m" };
        os.setf(os.fixed);
        os.precision(20);
        static_assert(n_samples > 0);
        println(os, '{');
        for (unsigned long i = 0; i < particles.size() - 1; ++i) {
            println(os, "    ", particles[i], ", ");
        }
        println(os, "    ", particles.back());
        println(os, '}');
        os.close();
    }
}
#if 0
TEST_CASE("Test LibPIC::RelativisticVDF::MaxwellianVDF::delta_f", "[LibPIC::RelativisticVDF::MaxwellianVDF::delta_f]")
{
    Real const O0 = 1., op = 4 * O0, c = op, beta1_eq = .1, T2OT1_eq = 5.35;
    Real const xi = .876, xiD1q1max = M_PI_2 * 0.8, psd_refresh_frequency = 0;
    long const q1min = -7, q1max = 15;
    auto const D1      = xiD1q1max / (xi * std::max(std::abs(q1min), std::abs(q1max)));
    auto const geo     = Geometry{ xi, D1, O0 };
    auto const kinetic = KineticPlasmaDesc{ { -O0, op }, 10, ShapeOrder::CIC, psd_refresh_frequency, ParticleScheme::delta_f, .001, 2.1 };
    auto const desc    = BiMaxPlasmaDesc(kinetic, beta1_eq, T2OT1_eq);
    auto const vdf     = RelativisticMaxwellianVDF(desc, geo, { q1min, q1max - q1min }, c);

    auto const g_desc = BiMaxPlasmaDesc({ { -O0, op }, 10, ShapeOrder::CIC }, beta1_eq * desc.marker_temp_ratio, T2OT1_eq);
    auto const g_vdf  = RelativisticMaxwellianVDF(g_desc, geo, { q1min, q1max - q1min }, c);

    CHECK(serialize(desc) == serialize(vdf.plasma_desc()));

    // check equilibrium macro variables
    CHECK(vdf.Nrefcell_div_Ntotal() == Approx{ 0.0796995979174554 }.epsilon(1e-10));

    std::array const etas{
        0.4287882239954878, 0.49761722933466956, 0.5815167468245395, 0.6800561885442317, 0.788001590917828,
        0.8920754659248894, 0.9704416860618577, 1., 0.9704416860618577, 0.8920754659248894, 0.788001590917828,
        0.6800561885442317, 0.5815167468245395, 0.49761722933466956, 0.4287882239954878, 0.3733660301900688,
        0.32911727304544663, 0.29391455144617107, 0.26596447719479277, 0.24383880547267717, 0.22643378958243887,
        0.21291237930064547, 0.2026501794964986
    };
    long q1 = q1min;
    for (auto const &eta : etas) {
        CurviCoord const pos{ Real(q1++) };

        auto const n0_ref = eta;
        auto const n0     = vdf.n0(pos);
        CHECK(*n0 == Approx{ n0_ref }.epsilon(1e-10));

        auto const nV0 = geo.cart_to_fac(vdf.nV0(pos), pos);
        CHECK(nV0.x == Approx{ 0 }.margin(1e-10));
        CHECK(nV0.y == Approx{ 0 }.margin(1e-10));
        CHECK(nV0.z == Approx{ 0 }.margin(1e-10));

        auto const nuv0 = geo.cart_to_fac(vdf.nuv0(pos), pos);
        REQUIRE(nuv0.tt == Approx{ *n0 * (c * c + .5 * beta1_eq * (.5 + T2OT1_eq * eta)) }.epsilon(1e-10));
        REQUIRE(nuv0.ts.x == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ts.y == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ts.z == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.xx == Approx{ *n0 * beta1_eq / 2 * (1 - beta1_eq / (c * c) * (3. / 4 + .5 * T2OT1_eq * eta)) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.yy == Approx{ *n0 * T2OT1_eq * eta * beta1_eq / 2 * (1 - beta1_eq / (c * c) * (1. / 4 + T2OT1_eq * eta)) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.zz == Approx{ *n0 * T2OT1_eq * eta * beta1_eq / 2 * (1 - beta1_eq / (c * c) * (1. / 4 + T2OT1_eq * eta)) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.xy == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.yz == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.zx == Approx{ 0 }.margin(1e-10));
    }

    // sampling
    auto const n_samples = 100000U;
    auto const particles = vdf.emit(n_samples);

    static_assert(n_samples > 100);
    std::for_each_n(begin(particles), 100, [&](Particle const &ptl) {
        REQUIRE(ptl.psd.weight == desc.initial_weight);
        REQUIRE(ptl.psd.marker == g_vdf.f0(ptl));
        REQUIRE(ptl.psd.real_f == Approx{ vdf.f0(ptl) + ptl.psd.weight * ptl.psd.marker }.epsilon(1e-10));
        REQUIRE(vdf.real_f0(ptl) == vdf.f0(ptl));
        REQUIRE(ptl.gamma == Approx{ std::sqrt(1 + dot(ptl.g_vel, ptl.g_vel) / (c * c)) }.epsilon(1e-10));

        auto const n     = *vdf.n0(ptl.pos);
        auto const vth1  = std::sqrt(beta1_eq);
        auto const vth2  = vth1 * std::sqrt(T2OT1_eq * n);
        auto const v1    = ptl.g_vel.x;
        auto const v2    = std::sqrt(ptl.g_vel.y * ptl.g_vel.y + ptl.g_vel.z * ptl.g_vel.z);
        auto const f_ref = n * std::exp(-v1 * v1 / (vth1 * vth1) - v2 * v2 / (vth2 * vth2))
                         / (4 * M_PI_2 / M_2_SQRTPI * vth1 * vth2 * vth2);
        REQUIRE(vdf.f0(ptl) == Approx{ f_ref }.epsilon(1e-10));
        auto const g_ref = n * std::exp(-v1 * v1 / (vth1 * vth1 * desc.marker_temp_ratio) - v2 * v2 / (vth2 * vth2 * desc.marker_temp_ratio))
                         / (4 * M_PI_2 / M_2_SQRTPI * vth1 * vth2 * vth2 * desc.marker_temp_ratio * std::sqrt(desc.marker_temp_ratio));
        REQUIRE(vdf.g0(ptl) == Approx{ g_ref }.epsilon(1e-10));
    });

    if constexpr (dump_samples) {
        std::ofstream os{ "/Users/kyungguk/Downloads/RelativisticMaxwellianVDF-delta_f.m" };
        os.setf(os.fixed);
        os.precision(20);
        static_assert(n_samples > 0);
        println(os, '{');
        for (unsigned long i = 0; i < particles.size() - 1; ++i) {
            println(os, "    ", particles[i], ", ");
        }
        println(os, "    ", particles.back());
        println(os, '}');
        os.close();
    }
}

TEST_CASE("Test LibPIC::RelativisticVDF::LossconeVDF::BiMax::Homogeneous", "[LibPIC::RelativisticVDF::LossconeVDF::BiMax::Homogeneous]")
{
    Real const O0 = 1., op = 4 * O0, c = op, beta1_eq = .1, T2OT1_eq = 5.35;
    Real const xi = 0, D1 = 1;
    long const q1min = -7, q1max = 15;
    auto const geo  = Geometry{ xi, D1, O0 };
    auto const desc = BiMaxPlasmaDesc({ { -O0, op }, 10, ShapeOrder::CIC }, beta1_eq, T2OT1_eq);
    auto const vdf  = RelativisticLossconeVDF(LossconePlasmaDesc{ desc }, geo, { q1min, q1max - q1min }, c);

    CHECK(serialize(LossconePlasmaDesc{ desc }) == serialize(vdf.plasma_desc()));

    // check equilibrium macro variables
    CHECK(vdf.Nrefcell_div_Ntotal() == Approx{ 1.0 / (q1max - q1min) }.epsilon(1e-10));

    for (long q1 = q1min; q1 <= q1max; ++q1) {
        CurviCoord const pos{ Real(q1) };
        auto const       eta               = 1;
        auto const       eta_b             = 1;
        auto const       beta_eq           = 1e-5;
        auto const       beta              = beta_eq * eta_b / eta;
        auto const       vth_ratio_squared = T2OT1_eq / (1 + beta_eq) * eta;

        auto const n0_ref = (eta - beta_eq * eta_b) / (1 - beta_eq);
        auto const n0     = vdf.n0(pos);
        REQUIRE(*n0 == Approx{ n0_ref }.epsilon(1e-10));

        auto const nV0 = geo.cart_to_fac(vdf.nV0(pos), pos);
        REQUIRE(nV0.x == Approx{ 0 }.margin(1e-10));
        REQUIRE(nV0.y == Approx{ 0 }.margin(1e-10));
        REQUIRE(nV0.z == Approx{ 0 }.margin(1e-10));

        auto const nuv0 = geo.cart_to_fac(vdf.nuv0(pos), pos);
        REQUIRE(nuv0.tt == Approx{ *n0 * (c * c + .5 * beta1_eq * (.5 + vth_ratio_squared * (1 + beta))) }.epsilon(1e-10));
        REQUIRE(nuv0.ts.x == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ts.y == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ts.z == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.xx == Approx{ *n0 * beta1_eq / 2 * (1 - beta1_eq / (c * c) * (3. / 4 + .5 * vth_ratio_squared * (1 + beta))) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.yy == Approx{ *n0 * vth_ratio_squared * (1 + beta) * beta1_eq / 2 * (1 - beta1_eq / (c * c) * (1. / 4 + vth_ratio_squared * (1 + (1 + beta) * beta) / (1 + beta))) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.zz == Approx{ *n0 * vth_ratio_squared * (1 + beta) * beta1_eq / 2 * (1 - beta1_eq / (c * c) * (1. / 4 + vth_ratio_squared * (1 + (1 + beta) * beta) / (1 + beta))) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.xy == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.yz == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.zx == Approx{ 0 }.margin(1e-10));
    }

    // sampling
    auto const n_samples = 100000U;
    auto const particles = vdf.emit(n_samples);

    static_assert(n_samples > 100);
    std::for_each_n(begin(particles), 100, [&](Particle const &ptl) {
        REQUIRE(ptl.psd.weight == 1);
        REQUIRE(ptl.psd.marker == vdf.f0(ptl));
        REQUIRE(ptl.psd.real_f == Approx{ ptl.psd.weight * ptl.psd.marker }.epsilon(1e-10));
        REQUIRE(vdf.real_f0(ptl) == vdf.f0(ptl));
        REQUIRE(ptl.gamma == Approx{ std::sqrt(1 + dot(ptl.g_vel, ptl.g_vel) / (c * c)) }.epsilon(1e-10));

        auto const B_div_B0          = geo.Bmag_div_B0(ptl.pos);
        auto const beta_eq           = 1e-5;
        auto const vth_ratio_squared = T2OT1_eq / (1 + beta_eq);
        auto const eta               = std::pow((1 - 1 / B_div_B0) * vth_ratio_squared + 1 / B_div_B0, -1);
        auto const eta_b             = std::pow((1 - 1 / B_div_B0) * beta_eq * vth_ratio_squared + 1 / B_div_B0, -1);
        auto const beta              = beta_eq * eta_b / eta;
        auto const vth1              = std::sqrt(beta1_eq);
        auto const vth2              = vth1 * std::sqrt(vth_ratio_squared * eta);
        auto const v1                = ptl.g_vel.x;
        auto const v2                = std::sqrt(ptl.g_vel.y * ptl.g_vel.y + ptl.g_vel.z * ptl.g_vel.z);
        auto const f_ref             = eta * (1 - beta) / (1 - beta_eq)
                         * std::exp(-v1 * v1 / (vth1 * vth1))
                         * (std::exp(-v2 * v2 / (vth2 * vth2)) - std::exp(-v2 * v2 / (beta * vth2 * vth2)))
                         / (4 * M_PI_2 / M_2_SQRTPI * vth1 * vth2 * vth2) / (1 - beta);
        REQUIRE(vdf.f0(ptl) == Approx{ f_ref }.epsilon(1e-10));
        auto const g_ref = (eta - beta_eq * eta_b) / (1 - beta_eq)
                         * std::exp(-v1 * v1 / (desc.marker_temp_ratio * vth1 * vth1))
                         * (std::exp(-v2 * v2 / (desc.marker_temp_ratio * vth2 * vth2)) - std::exp(-v2 * v2 / (beta * desc.marker_temp_ratio * vth2 * vth2)))
                         / (4 * M_PI_2 / M_2_SQRTPI * vth1 * vth2 * vth2 * desc.marker_temp_ratio * std::sqrt(desc.marker_temp_ratio)) / (1 - beta);
        REQUIRE(vdf.g0(ptl) == Approx{ g_ref }.epsilon(1e-10));
    });

    if constexpr (dump_samples) {
        std::ofstream os{ "/Users/kyungguk/Downloads/RelativisticLossconeVDF_BiMax-homogeneous.m" };
        os.setf(os.fixed);
        os.precision(20);
        static_assert(n_samples > 0);
        println(os, '{');
        for (unsigned long i = 0; i < particles.size() - 1; ++i) {
            println(os, "    ", particles[i], ", ");
        }
        println(os, "    ", particles.back());
        println(os, '}');
        os.close();
    }
}

TEST_CASE("Test LibPIC::RelativisticVDF::LossconeVDF::BiMax::Inhomogeneous", "[LibPIC::RelativisticVDF::LossconeVDF::BiMax::Inhomogeneous]")
{
    Real const O0 = 1., op = 4 * O0, c = op, beta1_eq = .1, T2OT1_eq = 5.35;
    Real const xi = .876, xiD1q1max = M_PI_2 * 0.8;
    long const q1min = -7, q1max = 15;
    auto const D1   = xiD1q1max / (xi * std::max(std::abs(q1min), std::abs(q1max)));
    auto const geo  = Geometry{ xi, D1, O0 };
    auto const desc = BiMaxPlasmaDesc({ { -O0, op }, 10, ShapeOrder::CIC }, beta1_eq, T2OT1_eq);
    auto const vdf  = RelativisticLossconeVDF(LossconePlasmaDesc{ desc }, geo, { q1min, q1max - q1min }, c);

    CHECK(serialize(LossconePlasmaDesc{ desc }) == serialize(vdf.plasma_desc()));

    // check equilibrium macro variables
    CHECK(vdf.Nrefcell_div_Ntotal() == Approx{ 0.07970134433386468 }.epsilon(1e-10));

    std::array const etas{
        0.42879123633025984, 0.49762030396573637, 0.5815197397970416, 0.6800588645090684, 0.7880036454881502, 0.8920766500089541,
        0.970442038846314, 1., 0.970442038846314, 0.8920766500089541, 0.7880036454881502, 0.6800588645090684, 0.5815197397970416,
        0.49762030396573637, 0.42879123633025984, 0.37336890766975, 0.3291199886157569, 0.29391710380836455, 0.2659668782647602,
        0.24384107315089862, 0.22643594386685753, 0.21291444034995088, 0.20265216678232262
    };
    std::array const eta_bs{
        1.4413912093052352, 1.3022091219893086, 1.1982158749947158, 1.1212616249880305, 1.0659200692069175, 1.0286058654348758,
        1.0070509750175554, 1., 1.0070509750175554, 1.0286058654348758, 1.0659200692069175, 1.1212616249880305, 1.1982158749947158,
        1.3022091219893086, 1.4413912093052352, 1.6281445589143488, 1.881749612870893, 2.2333127829938704, 2.7354240450511313,
        3.482518706776325, 4.657982242968423, 6.657178178048698, 10.466831843938639
    };
    long q1 = q1min;
    for (unsigned long i = 0; i < std::size(etas); ++i) {
        CurviCoord const pos{ Real(q1++) };

        auto const eta               = etas.at(i);
        auto const eta_b             = eta_bs.at(i);
        auto const beta_eq           = 1e-5;
        auto const beta              = beta_eq * eta_b / eta;
        auto const vth_ratio_squared = T2OT1_eq / (1 + beta_eq) * eta;

        auto const n0_ref = (eta - beta_eq * eta_b) / (1 - beta_eq);
        auto const n0     = vdf.n0(pos);
        REQUIRE(*n0 == Approx{ n0_ref }.epsilon(1e-10));

        auto const nV0 = geo.cart_to_fac(vdf.nV0(pos), pos);
        REQUIRE(nV0.x == Approx{ 0 }.margin(1e-10));
        REQUIRE(nV0.y == Approx{ 0 }.margin(1e-10));
        REQUIRE(nV0.z == Approx{ 0 }.margin(1e-10));

        auto const nuv0 = geo.cart_to_fac(vdf.nuv0(pos), pos);
        REQUIRE(nuv0.tt == Approx{ *n0 * (c * c + .5 * beta1_eq * (.5 + vth_ratio_squared * (1 + beta))) }.epsilon(1e-10));
        REQUIRE(nuv0.ts.x == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ts.y == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ts.z == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.xx == Approx{ *n0 * beta1_eq / 2 * (1 - beta1_eq / (c * c) * (3. / 4 + .5 * vth_ratio_squared * (1 + beta))) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.yy == Approx{ *n0 * vth_ratio_squared * (1 + beta) * beta1_eq / 2 * (1 - beta1_eq / (c * c) * (1. / 4 + vth_ratio_squared * (1 + (1 + beta) * beta) / (1 + beta))) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.zz == Approx{ *n0 * vth_ratio_squared * (1 + beta) * beta1_eq / 2 * (1 - beta1_eq / (c * c) * (1. / 4 + vth_ratio_squared * (1 + (1 + beta) * beta) / (1 + beta))) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.xy == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.yz == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.zx == Approx{ 0 }.margin(1e-10));
    }

    // sampling
    auto const n_samples = 100000U;
    auto const particles = vdf.emit(n_samples);

    static_assert(n_samples > 100);
    std::for_each_n(begin(particles), 100, [&](Particle const &ptl) {
        REQUIRE(ptl.psd.weight == 1);
        REQUIRE(ptl.psd.marker == vdf.f0(ptl));
        REQUIRE(ptl.psd.real_f == Approx{ ptl.psd.weight * ptl.psd.marker }.epsilon(1e-10));
        REQUIRE(vdf.real_f0(ptl) == vdf.f0(ptl));
        REQUIRE(ptl.gamma == Approx{ std::sqrt(1 + dot(ptl.g_vel, ptl.g_vel) / (c * c)) }.epsilon(1e-10));

        auto const B_div_B0          = geo.Bmag_div_B0(ptl.pos);
        auto const beta_eq           = 1e-5;
        auto const vth_ratio_squared = T2OT1_eq / (1 + beta_eq);
        auto const eta               = std::pow((1 - 1 / B_div_B0) * vth_ratio_squared + 1 / B_div_B0, -1);
        auto const eta_b             = std::pow((1 - 1 / B_div_B0) * beta_eq * vth_ratio_squared + 1 / B_div_B0, -1);
        auto const beta              = beta_eq * eta_b / eta;
        auto const vth1              = std::sqrt(beta1_eq);
        auto const vth2              = vth1 * std::sqrt(vth_ratio_squared * eta);
        auto const v1                = ptl.g_vel.x;
        auto const v2                = std::sqrt(ptl.g_vel.y * ptl.g_vel.y + ptl.g_vel.z * ptl.g_vel.z);
        auto const f_ref             = eta * (1 - beta) / (1 - beta_eq)
                         * std::exp(-v1 * v1 / (vth1 * vth1))
                         * (std::exp(-v2 * v2 / (vth2 * vth2)) - std::exp(-v2 * v2 / (beta * vth2 * vth2)))
                         / (4 * M_PI_2 / M_2_SQRTPI * vth1 * vth2 * vth2) / (1 - beta);
        REQUIRE(vdf.f0(ptl) == Approx{ f_ref }.epsilon(1e-10));
        auto const g_ref = (eta - beta_eq * eta_b) / (1 - beta_eq)
                         * std::exp(-v1 * v1 / (desc.marker_temp_ratio * vth1 * vth1))
                         * (std::exp(-v2 * v2 / (desc.marker_temp_ratio * vth2 * vth2)) - std::exp(-v2 * v2 / (beta * desc.marker_temp_ratio * vth2 * vth2)))
                         / (4 * M_PI_2 / M_2_SQRTPI * vth1 * vth2 * vth2 * desc.marker_temp_ratio * std::sqrt(desc.marker_temp_ratio)) / (1 - beta);
        REQUIRE(vdf.g0(ptl) == Approx{ g_ref }.epsilon(1e-10));
    });

    if constexpr (dump_samples) {
        std::ofstream os{ "/Users/kyungguk/Downloads/RelativisticLossconeVDF_BiMax-inhomogeneous.m" };
        os.setf(os.fixed);
        os.precision(20);
        static_assert(n_samples > 0);
        println(os, '{');
        for (unsigned long i = 0; i < particles.size() - 1; ++i) {
            println(os, "    ", particles[i], ", ");
        }
        println(os, "    ", particles.back());
        println(os, '}');
        os.close();
    }
}

TEST_CASE("Test LibPIC::RelativisticVDF::LossconeVDF::BiMax::delta_f", "[LibPIC::RelativisticVDF::LossconeVDF::BiMax::delta_f]")
{
    Real const O0 = 1., op = 4 * O0, c = op, beta1_eq = .1, T2OT1_eq = 5.35;
    Real const xi = .876, xiD1q1max = M_PI_2 * 0.8, psd_refresh_frequency = 0;
    long const q1min = -7, q1max = 15;
    auto const D1      = xiD1q1max / (xi * std::max(std::abs(q1min), std::abs(q1max)));
    auto const geo     = Geometry{ xi, D1, O0 };
    auto const kinetic = KineticPlasmaDesc{ { -O0, op }, 10, ShapeOrder::CIC, psd_refresh_frequency, ParticleScheme::delta_f, .001, 2.1 };
    auto const desc    = BiMaxPlasmaDesc(kinetic, beta1_eq, T2OT1_eq);
    auto const vdf     = RelativisticLossconeVDF(LossconePlasmaDesc{ desc }, geo, { q1min, q1max - q1min }, c);

    auto const g_desc = BiMaxPlasmaDesc({ { -O0, op }, 10, ShapeOrder::CIC }, beta1_eq * desc.marker_temp_ratio, T2OT1_eq);
    auto const g_vdf  = RelativisticLossconeVDF(LossconePlasmaDesc{ g_desc }, geo, { q1min, q1max - q1min }, c);

    CHECK(serialize(LossconePlasmaDesc{ desc }) == serialize(vdf.plasma_desc()));

    // check equilibrium macro variables
    CHECK(vdf.Nrefcell_div_Ntotal() == Approx{ 0.07970134433386468 }.epsilon(1e-10));

    std::array const etas{
        0.42879123633025984, 0.49762030396573637, 0.5815197397970416, 0.6800588645090684, 0.7880036454881502, 0.8920766500089541,
        0.970442038846314, 1., 0.970442038846314, 0.8920766500089541, 0.7880036454881502, 0.6800588645090684, 0.5815197397970416,
        0.49762030396573637, 0.42879123633025984, 0.37336890766975, 0.3291199886157569, 0.29391710380836455, 0.2659668782647602,
        0.24384107315089862, 0.22643594386685753, 0.21291444034995088, 0.20265216678232262
    };
    std::array const eta_bs{
        1.4413912093052352, 1.3022091219893086, 1.1982158749947158, 1.1212616249880305, 1.0659200692069175, 1.0286058654348758,
        1.0070509750175554, 1., 1.0070509750175554, 1.0286058654348758, 1.0659200692069175, 1.1212616249880305, 1.1982158749947158,
        1.3022091219893086, 1.4413912093052352, 1.6281445589143488, 1.881749612870893, 2.2333127829938704, 2.7354240450511313,
        3.482518706776325, 4.657982242968423, 6.657178178048698, 10.466831843938639
    };
    long q1 = q1min;
    for (unsigned long i = 0; i < std::size(etas); ++i) {
        CurviCoord const pos{ Real(q1++) };

        auto const eta               = etas.at(i);
        auto const eta_b             = eta_bs.at(i);
        auto const beta_eq           = 1e-5;
        auto const beta              = beta_eq * eta_b / eta;
        auto const vth_ratio_squared = T2OT1_eq / (1 + beta_eq) * eta;

        auto const n0_ref = (eta - beta_eq * eta_b) / (1 - beta_eq);
        auto const n0     = vdf.n0(pos);
        REQUIRE(*n0 == Approx{ n0_ref }.epsilon(1e-10));

        auto const nV0 = geo.cart_to_fac(vdf.nV0(pos), pos);
        REQUIRE(nV0.x == Approx{ 0 }.margin(1e-10));
        REQUIRE(nV0.y == Approx{ 0 }.margin(1e-10));
        REQUIRE(nV0.z == Approx{ 0 }.margin(1e-10));

        auto const nuv0 = geo.cart_to_fac(vdf.nuv0(pos), pos);
        REQUIRE(nuv0.tt == Approx{ *n0 * (c * c + .5 * beta1_eq * (.5 + vth_ratio_squared * (1 + beta))) }.epsilon(1e-10));
        REQUIRE(nuv0.ts.x == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ts.y == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ts.z == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.xx == Approx{ *n0 * beta1_eq / 2 * (1 - beta1_eq / (c * c) * (3. / 4 + .5 * vth_ratio_squared * (1 + beta))) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.yy == Approx{ *n0 * vth_ratio_squared * (1 + beta) * beta1_eq / 2 * (1 - beta1_eq / (c * c) * (1. / 4 + vth_ratio_squared * (1 + (1 + beta) * beta) / (1 + beta))) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.zz == Approx{ *n0 * vth_ratio_squared * (1 + beta) * beta1_eq / 2 * (1 - beta1_eq / (c * c) * (1. / 4 + vth_ratio_squared * (1 + (1 + beta) * beta) / (1 + beta))) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.xy == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.yz == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.zx == Approx{ 0 }.margin(1e-10));
    }

    // sampling
    auto const n_samples = 100000U;
    auto const particles = vdf.emit(n_samples);

    static_assert(n_samples > 100);
    std::for_each_n(begin(particles), 100, [&](Particle const &ptl) {
        REQUIRE(ptl.psd.weight == desc.initial_weight);
        REQUIRE(ptl.psd.marker == g_vdf.f0(ptl));
        REQUIRE(ptl.psd.real_f == Approx{ vdf.f0(ptl) + ptl.psd.weight * ptl.psd.marker }.epsilon(1e-10));
        REQUIRE(vdf.real_f0(ptl) == vdf.f0(ptl));
        REQUIRE(ptl.gamma == Approx{ std::sqrt(1 + dot(ptl.g_vel, ptl.g_vel) / (c * c)) }.epsilon(1e-10));

        auto const B_div_B0          = geo.Bmag_div_B0(ptl.pos);
        auto const beta_eq           = 1e-5;
        auto const vth_ratio_squared = T2OT1_eq / (1 + beta_eq);
        auto const eta               = std::pow((1 - 1 / B_div_B0) * vth_ratio_squared + 1 / B_div_B0, -1);
        auto const eta_b             = std::pow((1 - 1 / B_div_B0) * beta_eq * vth_ratio_squared + 1 / B_div_B0, -1);
        auto const beta              = beta_eq * eta_b / eta;
        auto const vth1              = std::sqrt(beta1_eq);
        auto const vth2              = vth1 * std::sqrt(vth_ratio_squared * eta);
        auto const v1                = ptl.g_vel.x;
        auto const v2                = std::sqrt(ptl.g_vel.y * ptl.g_vel.y + ptl.g_vel.z * ptl.g_vel.z);
        auto const f_ref             = eta * (1 - beta) / (1 - beta_eq)
                         * std::exp(-v1 * v1 / (vth1 * vth1))
                         * (std::exp(-v2 * v2 / (vth2 * vth2)) - std::exp(-v2 * v2 / (beta * vth2 * vth2)))
                         / (4 * M_PI_2 / M_2_SQRTPI * vth1 * vth2 * vth2) / (1 - beta);
        REQUIRE(vdf.f0(ptl) == Approx{ f_ref }.epsilon(1e-10));
        auto const g_ref = (eta - beta_eq * eta_b) / (1 - beta_eq)
                         * std::exp(-v1 * v1 / (desc.marker_temp_ratio * vth1 * vth1))
                         * (std::exp(-v2 * v2 / (desc.marker_temp_ratio * vth2 * vth2)) - std::exp(-v2 * v2 / (beta * desc.marker_temp_ratio * vth2 * vth2)))
                         / (4 * M_PI_2 / M_2_SQRTPI * vth1 * vth2 * vth2 * desc.marker_temp_ratio * std::sqrt(desc.marker_temp_ratio)) / (1 - beta);
        REQUIRE(vdf.g0(ptl) == Approx{ g_ref }.epsilon(1e-10));
    });

    if constexpr (dump_samples) {
        std::ofstream os{ "/Users/kyungguk/Downloads/RelativisticLossconeVDF_BiMax-delta_f.m" };
        os.setf(os.fixed);
        os.precision(20);
        static_assert(n_samples > 0);
        println(os, '{');
        for (unsigned long i = 0; i < particles.size() - 1; ++i) {
            println(os, "    ", particles[i], ", ");
        }
        println(os, "    ", particles.back());
        println(os, '}');
        os.close();
    }
}

TEST_CASE("Test LibPIC::RelativisticVDF::LossconeVDF::Loss::Homogeneous", "[LibPIC::RelativisticVDF::LossconeVDF::Loss::Homogeneous]")
{
    Real const O0 = 1., op = 4 * O0, c = op, beta1_eq = .1, T2OT1_eq = 5.35, beta_eq = .9;
    Real const xi = 0, D1 = 1;
    long const q1min = -7, q1max = 15;
    auto const geo     = Geometry{ xi, D1, O0 };
    auto const kinetic = KineticPlasmaDesc{ { -O0, op }, 10, ShapeOrder::CIC };
    auto const desc    = LossconePlasmaDesc(kinetic, beta1_eq, T2OT1_eq / (1 + beta_eq), beta_eq);
    auto const vdf     = RelativisticLossconeVDF(desc, geo, { q1min, q1max - q1min }, c);

    CHECK(serialize(desc) == serialize(vdf.plasma_desc()));

    // check equilibrium macro variables
    CHECK(vdf.Nrefcell_div_Ntotal() == Approx{ 1.0 / (q1max - q1min) }.epsilon(1e-10));

    for (long q1 = q1min; q1 <= q1max; ++q1) {
        CurviCoord const pos{ Real(q1) };
        auto const       eta               = 1;
        auto const       eta_b             = 1;
        auto const       beta              = beta_eq * eta_b / eta;
        auto const       vth_ratio_squared = T2OT1_eq / (1 + beta_eq) * eta;

        auto const n0_ref = (eta - beta_eq * eta_b) / (1 - beta_eq);
        auto const n0     = vdf.n0(pos);
        REQUIRE(*n0 == Approx{ n0_ref }.epsilon(1e-10));

        auto const nV0 = geo.cart_to_fac(vdf.nV0(pos), pos);
        REQUIRE(nV0.x == Approx{ 0 }.margin(1e-10));
        REQUIRE(nV0.y == Approx{ 0 }.margin(1e-10));
        REQUIRE(nV0.z == Approx{ 0 }.margin(1e-10));

        auto const nuv0 = geo.cart_to_fac(vdf.nuv0(pos), pos);
        REQUIRE(nuv0.tt == Approx{ *n0 * (c * c + .5 * beta1_eq * (.5 + vth_ratio_squared * (1 + beta))) }.epsilon(1e-10));
        REQUIRE(nuv0.ts.x == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ts.y == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ts.z == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.xx == Approx{ *n0 * beta1_eq / 2 * (1 - beta1_eq / (c * c) * (3. / 4 + .5 * vth_ratio_squared * (1 + beta))) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.yy == Approx{ *n0 * vth_ratio_squared * (1 + beta) * beta1_eq / 2 * (1 - beta1_eq / (c * c) * (1. / 4 + vth_ratio_squared * (1 + (1 + beta) * beta) / (1 + beta))) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.zz == Approx{ *n0 * vth_ratio_squared * (1 + beta) * beta1_eq / 2 * (1 - beta1_eq / (c * c) * (1. / 4 + vth_ratio_squared * (1 + (1 + beta) * beta) / (1 + beta))) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.xy == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.yz == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.zx == Approx{ 0 }.margin(1e-10));
    }

    // sampling
    auto const n_samples = 100000U;
    auto const particles = vdf.emit(n_samples);

    static_assert(n_samples > 100);
    std::for_each_n(begin(particles), 100, [&](Particle const &ptl) {
        REQUIRE(ptl.psd.weight == 1);
        REQUIRE(ptl.psd.marker == vdf.f0(ptl));
        REQUIRE(ptl.psd.real_f == Approx{ ptl.psd.weight * ptl.psd.marker }.epsilon(1e-10));
        REQUIRE(vdf.real_f0(ptl) == vdf.f0(ptl));
        REQUIRE(ptl.gamma == Approx{ std::sqrt(1 + dot(ptl.g_vel, ptl.g_vel) / (c * c)) }.epsilon(1e-10));

        auto const B_div_B0          = geo.Bmag_div_B0(ptl.pos);
        auto const vth_ratio_squared = T2OT1_eq / (1 + beta_eq);
        auto const eta               = std::pow((1 - 1 / B_div_B0) * vth_ratio_squared + 1 / B_div_B0, -1);
        auto const eta_b             = std::pow((1 - 1 / B_div_B0) * beta_eq * vth_ratio_squared + 1 / B_div_B0, -1);
        auto const beta              = beta_eq * eta_b / eta;
        auto const vth1              = std::sqrt(beta1_eq);
        auto const vth2              = vth1 * std::sqrt(vth_ratio_squared * eta);
        auto const v1                = ptl.g_vel.x;
        auto const v2                = std::sqrt(ptl.g_vel.y * ptl.g_vel.y + ptl.g_vel.z * ptl.g_vel.z);
        auto const f_ref             = eta * (1 - beta) / (1 - beta_eq)
                         * std::exp(-v1 * v1 / (vth1 * vth1))
                         * (std::exp(-v2 * v2 / (vth2 * vth2)) - std::exp(-v2 * v2 / (beta * vth2 * vth2)))
                         / (4 * M_PI_2 / M_2_SQRTPI * vth1 * vth2 * vth2) / (1 - beta);
        REQUIRE(vdf.f0(ptl) == Approx{ f_ref }.epsilon(1e-10));
        auto const g_ref = (eta - beta_eq * eta_b) / (1 - beta_eq)
                         * std::exp(-v1 * v1 / (desc.marker_temp_ratio * vth1 * vth1))
                         * (std::exp(-v2 * v2 / (desc.marker_temp_ratio * vth2 * vth2)) - std::exp(-v2 * v2 / (beta * desc.marker_temp_ratio * vth2 * vth2)))
                         / (4 * M_PI_2 / M_2_SQRTPI * vth1 * vth2 * vth2 * desc.marker_temp_ratio * std::sqrt(desc.marker_temp_ratio)) / (1 - beta);
        REQUIRE(vdf.g0(ptl) == Approx{ g_ref }.epsilon(1e-10));
    });

    if constexpr (dump_samples) {
        std::ofstream os{ "/Users/kyungguk/Downloads/RelativisticLossconeVDF_Loss-homogeneous.m" };
        os.setf(os.fixed);
        os.precision(20);
        static_assert(n_samples > 0);
        println(os, '{');
        for (unsigned long i = 0; i < particles.size() - 1; ++i) {
            println(os, "    ", particles[i], ", ");
        }
        println(os, "    ", particles.back());
        println(os, '}');
        os.close();
    }
}

TEST_CASE("Test LibPIC::RelativisticVDF::LossconeVDF::Loss::Inhomogeneous", "[LibPIC::RelativisticVDF::LossconeVDF::Loss::Inhomogeneous]")
{
    Real const O0 = 1., op = 4 * O0, c = op, beta1_eq = .1, T2OT1_eq = 5.35, beta_eq = .9;
    Real const xi = .876, xiD1q1max = M_PI_2 * 0.8;
    long const q1min = -7, q1max = 15;
    auto const D1      = xiD1q1max / (xi * std::max(std::abs(q1min), std::abs(q1max)));
    auto const geo     = Geometry{ xi, D1, O0 };
    auto const kinetic = KineticPlasmaDesc{ { -O0, op }, 10, ShapeOrder::CIC };
    auto const desc    = LossconePlasmaDesc(kinetic, beta1_eq, T2OT1_eq / (1 + beta_eq), beta_eq);
    auto const vdf     = RelativisticLossconeVDF(desc, geo, { q1min, q1max - q1min }, c);

    CHECK(serialize(desc) == serialize(vdf.plasma_desc()));

    // check equilibrium macro variables
    CHECK(vdf.Nrefcell_div_Ntotal() == Approx{ 0.09452655524436228 }.epsilon(1e-10));

    std::array const etas{
        0.64264404305085, 0.7035216274440845, 0.7689973601535578, 0.835852333403872, 0.8990376137203975, 0.9519272903927195,
        0.9874454952939664, 1., 0.9874454952939664, 0.9519272903927195, 0.8990376137203975, 0.835852333403872, 0.7689973601535578,
        0.7035216274440845, 0.64264404305085, 0.5880359901402682, 0.5402813355647949, 0.4993020809205341, 0.46467416286831936,
        0.4358328683746492, 0.4121939193295821, 0.39321849079530274, 0.37844426744095166
    };
    std::array const eta_bs{
        0.6803461521487998, 0.7374252176467319, 0.7975679579815613, 0.8576845064678151, 0.9133372658253722, 0.9590769728307393,
        0.9893716614043512, 1., 0.9893716614043512, 0.9590769728307393, 0.9133372658253722, 0.8576845064678151, 0.7975679579815613,
        0.7374252176467319, 0.6803461521487998, 0.6281659134692374, 0.5817544384045442, 0.5413336389508113, 0.5067409023782806,
        0.47761814536264413, 0.45353476372795026, 0.4340615580713739, 0.41881195104513863
    };
    long q1 = q1min;
    for (unsigned long i = 0; i < std::size(etas); ++i) {
        CurviCoord const pos{ Real(q1++) };

        auto const eta               = etas.at(i);
        auto const eta_b             = eta_bs.at(i);
        auto const beta              = beta_eq * eta_b / eta;
        auto const vth_ratio_squared = T2OT1_eq / (1 + beta_eq) * eta;

        auto const n0_ref = (eta - beta_eq * eta_b) / (1 - beta_eq);
        auto const n0     = vdf.n0(pos);
        REQUIRE(*n0 == Approx{ n0_ref }.epsilon(1e-10));

        auto const nV0 = geo.cart_to_fac(vdf.nV0(pos), pos);
        REQUIRE(nV0.x == Approx{ 0 }.margin(1e-10));
        REQUIRE(nV0.y == Approx{ 0 }.margin(1e-10));
        REQUIRE(nV0.z == Approx{ 0 }.margin(1e-10));

        auto const nuv0 = geo.cart_to_fac(vdf.nuv0(pos), pos);
        REQUIRE(nuv0.tt == Approx{ *n0 * (c * c + .5 * beta1_eq * (.5 + vth_ratio_squared * (1 + beta))) }.epsilon(1e-10));
        REQUIRE(nuv0.ts.x == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ts.y == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ts.z == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.xx == Approx{ *n0 * beta1_eq / 2 * (1 - beta1_eq / (c * c) * (3. / 4 + .5 * vth_ratio_squared * (1 + beta))) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.yy == Approx{ *n0 * vth_ratio_squared * (1 + beta) * beta1_eq / 2 * (1 - beta1_eq / (c * c) * (1. / 4 + vth_ratio_squared * (1 + (1 + beta) * beta) / (1 + beta))) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.zz == Approx{ *n0 * vth_ratio_squared * (1 + beta) * beta1_eq / 2 * (1 - beta1_eq / (c * c) * (1. / 4 + vth_ratio_squared * (1 + (1 + beta) * beta) / (1 + beta))) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.xy == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.yz == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.zx == Approx{ 0 }.margin(1e-10));
    }

    // sampling
    auto const n_samples = 100000U;
    auto const particles = vdf.emit(n_samples);

    static_assert(n_samples > 100);
    std::for_each_n(begin(particles), 100, [&](Particle const &ptl) {
        REQUIRE(ptl.psd.weight == 1);
        REQUIRE(ptl.psd.marker == vdf.f0(ptl));
        REQUIRE(ptl.psd.real_f == Approx{ ptl.psd.weight * ptl.psd.marker }.epsilon(1e-10));
        REQUIRE(vdf.real_f0(ptl) == vdf.f0(ptl));
        REQUIRE(ptl.gamma == Approx{ std::sqrt(1 + dot(ptl.g_vel, ptl.g_vel) / (c * c)) }.epsilon(1e-10));

        auto const B_div_B0          = geo.Bmag_div_B0(ptl.pos);
        auto const vth_ratio_squared = T2OT1_eq / (1 + beta_eq);
        auto const eta               = std::pow((1 - 1 / B_div_B0) * vth_ratio_squared + 1 / B_div_B0, -1);
        auto const eta_b             = std::pow((1 - 1 / B_div_B0) * beta_eq * vth_ratio_squared + 1 / B_div_B0, -1);
        auto const beta              = beta_eq * eta_b / eta;
        auto const vth1              = std::sqrt(beta1_eq);
        auto const vth2              = vth1 * std::sqrt(vth_ratio_squared * eta);
        auto const v1                = ptl.g_vel.x;
        auto const v2                = std::sqrt(ptl.g_vel.y * ptl.g_vel.y + ptl.g_vel.z * ptl.g_vel.z);
        auto const f_ref             = eta * (1 - beta) / (1 - beta_eq)
                         * std::exp(-v1 * v1 / (vth1 * vth1))
                         * (std::exp(-v2 * v2 / (vth2 * vth2)) - std::exp(-v2 * v2 / (beta * vth2 * vth2)))
                         / (4 * M_PI_2 / M_2_SQRTPI * vth1 * vth2 * vth2) / (1 - beta);
        REQUIRE(vdf.f0(ptl) == Approx{ f_ref }.epsilon(1e-10));
        auto const g_ref = (eta - beta_eq * eta_b) / (1 - beta_eq)
                         * std::exp(-v1 * v1 / (desc.marker_temp_ratio * vth1 * vth1))
                         * (std::exp(-v2 * v2 / (desc.marker_temp_ratio * vth2 * vth2)) - std::exp(-v2 * v2 / (beta * desc.marker_temp_ratio * vth2 * vth2)))
                         / (4 * M_PI_2 / M_2_SQRTPI * vth1 * vth2 * vth2 * desc.marker_temp_ratio * std::sqrt(desc.marker_temp_ratio)) / (1 - beta);
        REQUIRE(vdf.g0(ptl) == Approx{ g_ref }.epsilon(1e-10));
    });

    if constexpr (dump_samples) {
        std::ofstream os{ "/Users/kyungguk/Downloads/RelativisticLossconeVDF_Loss-inhomogeneous.m" };
        os.setf(os.fixed);
        os.precision(20);
        static_assert(n_samples > 0);
        println(os, '{');
        for (unsigned long i = 0; i < particles.size() - 1; ++i) {
            println(os, "    ", particles[i], ", ");
        }
        println(os, "    ", particles.back());
        println(os, '}');
        os.close();
    }
}

TEST_CASE("Test LibPIC::RelativisticVDF::LossconeVDF::Loss::delta_f", "[LibPIC::RelativisticVDF::LossconeVDF::Loss::delta_f]")
{
    Real const O0 = 1., op = 4 * O0, c = op, beta1_eq = .1, T2OT1_eq = 5.35, beta_eq = .9;
    Real const xi = .876, xiD1q1max = M_PI_2 * 0.8, psd_refresh_frequency = 0;
    long const q1min = -7, q1max = 15;
    auto const D1      = xiD1q1max / (xi * std::max(std::abs(q1min), std::abs(q1max)));
    auto const geo     = Geometry{ xi, D1, O0 };
    auto const kinetic = KineticPlasmaDesc{ { -O0, op }, 10, ShapeOrder::CIC, psd_refresh_frequency, ParticleScheme::delta_f, .001, 2.1 };
    auto const desc    = LossconePlasmaDesc(kinetic, beta1_eq, T2OT1_eq / (1 + beta_eq), beta_eq);
    auto const vdf     = RelativisticLossconeVDF(desc, geo, { q1min, q1max - q1min }, c);

    auto const g_desc = BiMaxPlasmaDesc({ { -O0, op }, 10, ShapeOrder::CIC }, beta1_eq * desc.marker_temp_ratio, T2OT1_eq);
    auto const g_vdf  = RelativisticLossconeVDF(LossconePlasmaDesc{ g_desc, beta_eq }, geo, { q1min, q1max - q1min }, c);

    CHECK(serialize(desc) == serialize(vdf.plasma_desc()));

    // check equilibrium macro variables
    CHECK(vdf.Nrefcell_div_Ntotal() == Approx{ 0.09452655524436228 }.epsilon(1e-10));

    std::array const etas{
        0.64264404305085, 0.7035216274440845, 0.7689973601535578, 0.835852333403872, 0.8990376137203975, 0.9519272903927195,
        0.9874454952939664, 1., 0.9874454952939664, 0.9519272903927195, 0.8990376137203975, 0.835852333403872, 0.7689973601535578,
        0.7035216274440845, 0.64264404305085, 0.5880359901402682, 0.5402813355647949, 0.4993020809205341, 0.46467416286831936,
        0.4358328683746492, 0.4121939193295821, 0.39321849079530274, 0.37844426744095166
    };
    std::array const eta_bs{
        0.6803461521487998, 0.7374252176467319, 0.7975679579815613, 0.8576845064678151, 0.9133372658253722, 0.9590769728307393,
        0.9893716614043512, 1., 0.9893716614043512, 0.9590769728307393, 0.9133372658253722, 0.8576845064678151, 0.7975679579815613,
        0.7374252176467319, 0.6803461521487998, 0.6281659134692374, 0.5817544384045442, 0.5413336389508113, 0.5067409023782806,
        0.47761814536264413, 0.45353476372795026, 0.4340615580713739, 0.41881195104513863
    };
    long q1 = q1min;
    for (unsigned long i = 0; i < std::size(etas); ++i) {
        CurviCoord const pos{ Real(q1++) };

        auto const eta               = etas.at(i);
        auto const eta_b             = eta_bs.at(i);
        auto const beta              = beta_eq * eta_b / eta;
        auto const vth_ratio_squared = T2OT1_eq / (1 + beta_eq) * eta;

        auto const n0_ref = (eta - beta_eq * eta_b) / (1 - beta_eq);
        auto const n0     = vdf.n0(pos);
        REQUIRE(*n0 == Approx{ n0_ref }.epsilon(1e-10));

        auto const nV0 = geo.cart_to_fac(vdf.nV0(pos), pos);
        REQUIRE(nV0.x == Approx{ 0 }.margin(1e-10));
        REQUIRE(nV0.y == Approx{ 0 }.margin(1e-10));
        REQUIRE(nV0.z == Approx{ 0 }.margin(1e-10));

        auto const nuv0 = geo.cart_to_fac(vdf.nuv0(pos), pos);
        REQUIRE(nuv0.tt == Approx{ *n0 * (c * c + .5 * beta1_eq * (.5 + vth_ratio_squared * (1 + beta))) }.epsilon(1e-10));
        REQUIRE(nuv0.ts.x == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ts.y == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ts.z == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.xx == Approx{ *n0 * beta1_eq / 2 * (1 - beta1_eq / (c * c) * (3. / 4 + .5 * vth_ratio_squared * (1 + beta))) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.yy == Approx{ *n0 * vth_ratio_squared * (1 + beta) * beta1_eq / 2 * (1 - beta1_eq / (c * c) * (1. / 4 + vth_ratio_squared * (1 + (1 + beta) * beta) / (1 + beta))) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.zz == Approx{ *n0 * vth_ratio_squared * (1 + beta) * beta1_eq / 2 * (1 - beta1_eq / (c * c) * (1. / 4 + vth_ratio_squared * (1 + (1 + beta) * beta) / (1 + beta))) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.xy == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.yz == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.zx == Approx{ 0 }.margin(1e-10));
    }

    // sampling
    auto const n_samples = 100000U;
    auto const particles = vdf.emit(n_samples);

    static_assert(n_samples > 100);
    std::for_each_n(begin(particles), 100, [&](Particle const &ptl) {
        REQUIRE(ptl.psd.weight == desc.initial_weight);
        REQUIRE(ptl.psd.marker == g_vdf.f0(ptl));
        REQUIRE(ptl.psd.real_f == Approx{ vdf.f0(ptl) + ptl.psd.weight * ptl.psd.marker }.epsilon(1e-10));
        REQUIRE(vdf.real_f0(ptl) == vdf.f0(ptl));
        REQUIRE(ptl.gamma == Approx{ std::sqrt(1 + dot(ptl.g_vel, ptl.g_vel) / (c * c)) }.epsilon(1e-10));

        auto const B_div_B0          = geo.Bmag_div_B0(ptl.pos);
        auto const vth_ratio_squared = T2OT1_eq / (1 + beta_eq);
        auto const eta               = std::pow((1 - 1 / B_div_B0) * vth_ratio_squared + 1 / B_div_B0, -1);
        auto const eta_b             = std::pow((1 - 1 / B_div_B0) * beta_eq * vth_ratio_squared + 1 / B_div_B0, -1);
        auto const beta              = beta_eq * eta_b / eta;
        auto const vth1              = std::sqrt(beta1_eq);
        auto const vth2              = vth1 * std::sqrt(vth_ratio_squared * eta);
        auto const v1                = ptl.g_vel.x;
        auto const v2                = std::sqrt(ptl.g_vel.y * ptl.g_vel.y + ptl.g_vel.z * ptl.g_vel.z);
        auto const f_ref             = eta * (1 - beta) / (1 - beta_eq)
                         * std::exp(-v1 * v1 / (vth1 * vth1))
                         * (std::exp(-v2 * v2 / (vth2 * vth2)) - std::exp(-v2 * v2 / (beta * vth2 * vth2)))
                         / (4 * M_PI_2 / M_2_SQRTPI * vth1 * vth2 * vth2) / (1 - beta);
        REQUIRE(vdf.f0(ptl) == Approx{ f_ref }.epsilon(1e-10));
        auto const g_ref = (eta - beta_eq * eta_b) / (1 - beta_eq)
                         * std::exp(-v1 * v1 / (desc.marker_temp_ratio * vth1 * vth1))
                         * (std::exp(-v2 * v2 / (desc.marker_temp_ratio * vth2 * vth2)) - std::exp(-v2 * v2 / (beta * desc.marker_temp_ratio * vth2 * vth2)))
                         / (4 * M_PI_2 / M_2_SQRTPI * vth1 * vth2 * vth2 * desc.marker_temp_ratio * std::sqrt(desc.marker_temp_ratio)) / (1 - beta);
        REQUIRE(vdf.g0(ptl) == Approx{ g_ref }.epsilon(1e-10));
    });

    if constexpr (dump_samples) {
        std::ofstream os{ "/Users/kyungguk/Downloads/RelativisticLossconeVDF_Loss-delta_f.m" };
        os.setf(os.fixed);
        os.precision(20);
        static_assert(n_samples > 0);
        println(os, '{');
        for (unsigned long i = 0; i < particles.size() - 1; ++i) {
            println(os, "    ", particles[i], ", ");
        }
        println(os, "    ", particles.back());
        println(os, '}');
        os.close();
    }
}

TEST_CASE("Test LibPIC::RelativisticVDF::PartialShellVDF::Homogeneous::Maxwellian", "[LibPIC::RelativisticVDF::PartialShellVDF::Homogeneous::Maxwellian]")
{
    Real const O0 = 1, op = 4 * O0, c = op, beta = 0.1;
    Real const xi = 0, D1 = 1;
    long const q1min = -7, q1max = 15;
    auto const geo  = Geometry{ xi, D1, O0 };
    auto const desc = PartialShellPlasmaDesc({ { -O0, op }, 10, ShapeOrder::CIC }, beta);
    auto const vdf  = RelativisticPartialShellVDF(desc, geo, { q1min, q1max - q1min }, c);

    CHECK(serialize(desc) == serialize(vdf.plasma_desc()));

    // check equilibrium macro variables
    CHECK(vdf.Nrefcell_div_Ntotal() == Approx{ 1.0 / (q1max - q1min) }.epsilon(1e-10));

    for (long q1 = q1min; q1 <= q1max; ++q1) {
        CurviCoord const pos{ Real(q1) };
        auto const       eta = 1;

        auto const n0_ref = eta;
        auto const n0     = vdf.n0(pos);
        CHECK(*n0 == Approx{ n0_ref }.epsilon(1e-10));

        auto const nV0 = geo.cart_to_fac(vdf.nV0(pos), pos);
        CHECK(nV0.x == Approx{ 0 }.margin(1e-10));
        CHECK(nV0.y == Approx{ 0 }.margin(1e-10));
        CHECK(nV0.z == Approx{ 0 }.margin(1e-10));

        auto const nuv0 = geo.cart_to_fac(vdf.nuv0(pos), pos);
        auto const xs   = desc.vs / std::sqrt(beta);
        auto const xs2  = xs * xs;
        auto const Ab   = .5 * (xs * std::exp(-xs2) + 2 / M_2_SQRTPI * (.5 + xs2) * std::erfc(-xs));
        auto const T    = .5 * desc.beta / Ab * (xs * (2.5 + xs2) * std::exp(-xs2) + 2 / M_2_SQRTPI * (.75 + xs2 * (3 + xs2)) * std::erfc(-xs));
        auto const TT   = .5 * desc.beta * desc.beta / Ab * (xs * (33. / 4 + xs2 * (7 + xs2)) * std::exp(-xs2) + 2 / M_2_SQRTPI * (15. / 8 + xs2 * (45. / 4 + 15. / 2 * xs2 + xs2 * xs2)) * std::erfc(-xs));
        REQUIRE(nuv0.tt == Approx{ *n0 * (c * c + .5 * T) }.epsilon(1e-10));
        REQUIRE(nuv0.ts.x == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ts.y == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ts.z == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.xx == Approx{ *n0 / Real(3 + desc.zeta) * (T - .5 * TT / (c * c)) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.yy == Approx{ *n0 * (1 + .5 * desc.zeta) / Real(3 + desc.zeta) * (T - .5 * TT / (c * c)) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.zz == Approx{ *n0 * (1 + .5 * desc.zeta) / Real(3 + desc.zeta) * (T - .5 * TT / (c * c)) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.xy == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.yz == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.zx == Approx{ 0 }.margin(1e-10));
    }

    // sampling
    auto const n_samples = 100000U;
    auto const particles = vdf.emit(n_samples);

    static_assert(n_samples > 100);
    std::for_each_n(begin(particles), 100, [&](Particle const &ptl) {
        REQUIRE(ptl.psd.weight == 1);
        REQUIRE(ptl.psd.marker == vdf.f0(ptl));
        REQUIRE(ptl.psd.real_f == Approx{ ptl.psd.weight * ptl.psd.marker }.epsilon(1e-10));
        REQUIRE(vdf.real_f0(ptl) == vdf.f0(ptl));
        REQUIRE(ptl.gamma == Approx{ std::sqrt(1 + dot(ptl.g_vel, ptl.g_vel) / (c * c)) }.epsilon(1e-10));

        auto const n     = *vdf.n0(ptl.pos);
        auto const vth2  = beta;
        auto const v     = std::sqrt(dot(ptl.g_vel, ptl.g_vel));
        auto const f_ref = n * std::exp(-v * v / vth2) / (M_PI * 2 / M_2_SQRTPI * vth2 * std::sqrt(vth2));
        REQUIRE(vdf.f0(ptl) == Approx{ f_ref }.epsilon(1e-10));

        auto const marker_vth2 = vth2 * desc.marker_temp_ratio;
        auto const g_ref       = n * std::exp(-v * v / marker_vth2) / (M_PI * 2 / M_2_SQRTPI * marker_vth2 * std::sqrt(marker_vth2));
        REQUIRE(vdf.g0(ptl) == Approx{ g_ref }.epsilon(1e-10));
    });

    if constexpr (dump_samples) {
        std::ofstream os{ "/Users/kyungguk/Downloads/RelativisticPartialShellVDF-homogeneous-maxwellian.m" };
        os.setf(os.fixed);
        os.precision(20);
        static_assert(n_samples > 0);
        println(os, '{');
        for (unsigned long i = 0; i < particles.size() - 1; ++i) {
            println(os, "    ", particles[i], ", ");
        }
        println(os, "    ", particles.back());
        println(os, '}');
        os.close();
    }
}

TEST_CASE("Test LibPIC::RelativisticVDF::PartialShellVDF::Homogeneous::IsotropicShell", "[LibPIC::RelativisticVDF::PartialShellVDF::Homogeneous::IsotropicShell]")
{
    Real const O0 = 1, op = 4 * O0, c = op, beta = 0.1, vs = 10;
    Real const xi = 0, D1 = 1, psd_refresh_frequency = 0;
    long const q1min = -7, q1max = 15;
    auto const geo  = Geometry{ xi, D1, O0 };
    auto const desc = PartialShellPlasmaDesc({ { -O0, op }, 10, ShapeOrder::CIC, psd_refresh_frequency, delta_f, .001, 2.1 }, beta, 0U, vs);
    auto const vdf  = RelativisticPartialShellVDF(desc, geo, { q1min, q1max - q1min }, c);

    auto const g_desc = PartialShellPlasmaDesc({ { -O0, op }, 10, ShapeOrder::CIC }, beta * desc.marker_temp_ratio, 0U, vs);
    auto const g_vdf  = RelativisticPartialShellVDF(g_desc, geo, { q1min, q1max - q1min }, c);

    CHECK(serialize(desc) == serialize(vdf.plasma_desc()));

    // check equilibrium macro variables
    CHECK(vdf.Nrefcell_div_Ntotal() == Approx{ 1.0 / (q1max - q1min) }.epsilon(1e-10));

    for (long q1 = q1min; q1 <= q1max; ++q1) {
        CurviCoord const pos{ Real(q1) };
        auto const       eta = 1;

        auto const n0_ref = eta;
        auto const n0     = vdf.n0(pos);
        CHECK(*n0 == Approx{ n0_ref }.epsilon(1e-10));

        auto const nV0 = geo.cart_to_fac(vdf.nV0(pos), pos);
        CHECK(nV0.x == Approx{ 0 }.margin(1e-10));
        CHECK(nV0.y == Approx{ 0 }.margin(1e-10));
        CHECK(nV0.z == Approx{ 0 }.margin(1e-10));

        auto const nuv0 = geo.cart_to_fac(vdf.nuv0(pos), pos);
        auto const xs   = desc.vs / std::sqrt(beta);
        auto const xs2  = xs * xs;
        auto const Ab   = .5 * (xs * std::exp(-xs2) + 2 / M_2_SQRTPI * (.5 + xs2) * std::erfc(-xs));
        auto const T    = .5 * desc.beta / Ab * (xs * (2.5 + xs2) * std::exp(-xs2) + 2 / M_2_SQRTPI * (.75 + xs2 * (3 + xs2)) * std::erfc(-xs));
        auto const TT   = .5 * desc.beta * desc.beta / Ab * (xs * (33. / 4 + xs2 * (7 + xs2)) * std::exp(-xs2) + 2 / M_2_SQRTPI * (15. / 8 + xs2 * (45. / 4 + 15. / 2 * xs2 + xs2 * xs2)) * std::erfc(-xs));
        REQUIRE(nuv0.tt == Approx{ *n0 * (c * c + .5 * T) }.epsilon(1e-10));
        REQUIRE(nuv0.ts.x == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ts.y == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ts.z == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.xx == Approx{ *n0 / Real(3 + desc.zeta) * (T - .5 * TT / (c * c)) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.yy == Approx{ *n0 * (1 + .5 * desc.zeta) / Real(3 + desc.zeta) * (T - .5 * TT / (c * c)) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.zz == Approx{ *n0 * (1 + .5 * desc.zeta) / Real(3 + desc.zeta) * (T - .5 * TT / (c * c)) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.xy == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.yz == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.zx == Approx{ 0 }.margin(1e-10));
    }

    // sampling
    auto const n_samples = 100000U;
    auto const particles = vdf.emit(n_samples);

    static_assert(n_samples > 100);
    std::for_each_n(begin(particles), 100, [&](Particle const &ptl) {
        REQUIRE(ptl.psd.weight == desc.initial_weight);
        REQUIRE(ptl.psd.marker == g_vdf.f0(ptl));
        REQUIRE(ptl.psd.real_f == Approx{ vdf.f0(ptl) + ptl.psd.weight * ptl.psd.marker }.epsilon(1e-10));
        REQUIRE(vdf.real_f0(ptl) == vdf.f0(ptl));
        REQUIRE(ptl.gamma == Approx{ std::sqrt(1 + dot(ptl.g_vel, ptl.g_vel) / (c * c)) }.epsilon(1e-10));

        auto const n = *vdf.n0(ptl.pos);
        {
            auto const vth2  = beta;
            auto const xs    = desc.vs / std::sqrt(vth2);
            auto const Ab    = .5 * (xs * std::exp(-xs * xs) + 2 / M_2_SQRTPI * (.5 + xs * xs) * std::erfc(-xs));
            auto const Bz    = 2 / M_2_SQRTPI * std::tgamma(1 + .5 * desc.zeta) / std::tgamma(1.5 + .5 * desc.zeta);
            auto const v     = std::sqrt(dot(ptl.g_vel, ptl.g_vel)) - desc.vs;
            auto const f_ref = n * std::exp(-v * v / vth2) / (M_PI * 2 * vth2 * std::sqrt(vth2) * Ab * Bz);
            REQUIRE(vdf.f0(ptl) == Approx{ f_ref }.epsilon(1e-10));
        }
        {
            auto const vth2  = beta * desc.marker_temp_ratio;
            auto const xs    = desc.vs / std::sqrt(vth2);
            auto const Ab    = .5 * (xs * std::exp(-xs * xs) + 2 / M_2_SQRTPI * (.5 + xs * xs) * std::erfc(-xs));
            auto const Bz    = 2 / M_2_SQRTPI * std::tgamma(1 + .5 * desc.zeta) / std::tgamma(1.5 + .5 * desc.zeta);
            auto const v     = std::sqrt(dot(ptl.g_vel, ptl.g_vel)) - desc.vs;
            auto const g_ref = n * std::exp(-v * v / vth2) / (M_PI * 2 * vth2 * std::sqrt(vth2) * Ab * Bz);
            REQUIRE(vdf.g0(ptl) == Approx{ g_ref }.epsilon(1e-10));
        }
    });

    if constexpr (dump_samples) {
        std::ofstream os{ "/Users/kyungguk/Downloads/RelativisticPartialShellVDF-homogeneous-isotropic_shell.m" };
        os.setf(os.fixed);
        os.precision(20);
        static_assert(n_samples > 0);
        println(os, '{');
        for (unsigned long i = 0; i < particles.size() - 1; ++i) {
            println(os, "    ", particles[i], ", ");
        }
        println(os, "    ", particles.back());
        println(os, '}');
        os.close();
    }
}

TEST_CASE("Test LibPIC::RelativisticVDF::PartialShellVDF::Homogeneous::AnisotropicShell", "[LibPIC::RelativisticVDF::PartialShellVDF::Homogeneous::AnisotropicShell]")
{
    unsigned const zeta = 30;
    Real const     O0 = 1, op = 4 * O0, c = op, beta = 0.1, vs = 10;
    Real const     xi = 0, D1 = 1, psd_refresh_frequency = 0;
    long const     q1min = -7, q1max = 15;
    auto const     geo  = Geometry{ xi, D1, O0 };
    auto const     desc = PartialShellPlasmaDesc({ { -O0, op }, 10, ShapeOrder::CIC, psd_refresh_frequency, delta_f, .001, 2.1 }, beta, zeta, vs);
    auto const     vdf  = RelativisticPartialShellVDF(desc, geo, { q1min, q1max - q1min }, c);

    auto const g_desc = PartialShellPlasmaDesc({ { -O0, op }, 10, ShapeOrder::CIC }, beta * desc.marker_temp_ratio, zeta, vs);
    auto const g_vdf  = RelativisticPartialShellVDF(g_desc, geo, { q1min, q1max - q1min }, c);

    CHECK(serialize(desc) == serialize(vdf.plasma_desc()));

    // check equilibrium macro variables
    CHECK(vdf.Nrefcell_div_Ntotal() == Approx{ 1.0 / (q1max - q1min) }.epsilon(1e-10));

    for (long q1 = q1min; q1 <= q1max; ++q1) {
        CurviCoord const pos{ Real(q1) };
        auto const       eta = 1;

        auto const n0_ref = eta;
        auto const n0     = vdf.n0(pos);
        CHECK(*n0 == Approx{ n0_ref }.epsilon(1e-10));

        auto const nV0 = geo.cart_to_fac(vdf.nV0(pos), pos);
        CHECK(nV0.x == Approx{ 0 }.margin(1e-10));
        CHECK(nV0.y == Approx{ 0 }.margin(1e-10));
        CHECK(nV0.z == Approx{ 0 }.margin(1e-10));

        auto const nuv0 = geo.cart_to_fac(vdf.nuv0(pos), pos);
        auto const xs   = desc.vs / std::sqrt(beta);
        auto const xs2  = xs * xs;
        auto const Ab   = .5 * (xs * std::exp(-xs2) + 2 / M_2_SQRTPI * (.5 + xs2) * std::erfc(-xs));
        auto const T    = .5 * desc.beta / Ab * (xs * (2.5 + xs2) * std::exp(-xs2) + 2 / M_2_SQRTPI * (.75 + xs2 * (3 + xs2)) * std::erfc(-xs));
        auto const TT   = .5 * desc.beta * desc.beta / Ab * (xs * (33. / 4 + xs2 * (7 + xs2)) * std::exp(-xs2) + 2 / M_2_SQRTPI * (15. / 8 + xs2 * (45. / 4 + 15. / 2 * xs2 + xs2 * xs2)) * std::erfc(-xs));
        REQUIRE(nuv0.tt == Approx{ *n0 * (c * c + .5 * T) }.epsilon(1e-10));
        REQUIRE(nuv0.ts.x == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ts.y == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ts.z == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.xx == Approx{ *n0 / Real(3 + desc.zeta) * (T - .5 * TT / (c * c)) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.yy == Approx{ *n0 * (1 + .5 * desc.zeta) / Real(3 + desc.zeta) * (T - .5 * TT / (c * c)) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.zz == Approx{ *n0 * (1 + .5 * desc.zeta) / Real(3 + desc.zeta) * (T - .5 * TT / (c * c)) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.xy == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.yz == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.zx == Approx{ 0 }.margin(1e-10));
    }

    // sampling
    auto const n_samples = 100000U;
    auto const particles = vdf.emit(n_samples);

    static_assert(n_samples > 100);
    std::for_each_n(begin(particles), 100, [&](Particle const &ptl) {
        REQUIRE(ptl.psd.weight == desc.initial_weight);
        REQUIRE(ptl.psd.marker == g_vdf.f0(ptl));
        REQUIRE(ptl.psd.real_f == Approx{ vdf.f0(ptl) + ptl.psd.weight * ptl.psd.marker }.epsilon(1e-10));
        REQUIRE(vdf.real_f0(ptl) == vdf.f0(ptl));
        REQUIRE(ptl.gamma == Approx{ std::sqrt(1 + dot(ptl.g_vel, ptl.g_vel) / (c * c)) }.epsilon(1e-10));

        auto const n = *vdf.n0(ptl.pos);
        {
            auto const vth2  = beta;
            auto const xs    = desc.vs / std::sqrt(vth2);
            auto const Ab    = .5 * (xs * std::exp(-xs * xs) + 2 / M_2_SQRTPI * (.5 + xs * xs) * std::erfc(-xs));
            auto const Bz    = 2 / M_2_SQRTPI * std::tgamma(1 + .5 * desc.zeta) / std::tgamma(1.5 + .5 * desc.zeta);
            auto const v     = std::sqrt(dot(ptl.g_vel, ptl.g_vel)) - desc.vs;
            auto const alpha = std::acos(ptl.g_vel.x / (v + desc.vs));
            auto const f_ref = n * std::exp(-v * v / vth2) * std::pow(std::sin(alpha), zeta) / (M_PI * 2 * vth2 * std::sqrt(vth2) * Ab * Bz);
            REQUIRE(vdf.f0(ptl) == Approx{ f_ref }.epsilon(1e-10));
        }
        {
            auto const vth2  = beta * desc.marker_temp_ratio;
            auto const xs    = desc.vs / std::sqrt(vth2);
            auto const Ab    = .5 * (xs * std::exp(-xs * xs) + 2 / M_2_SQRTPI * (.5 + xs * xs) * std::erfc(-xs));
            auto const Bz    = 2 / M_2_SQRTPI * std::tgamma(1 + .5 * desc.zeta) / std::tgamma(1.5 + .5 * desc.zeta);
            auto const v     = std::sqrt(dot(ptl.g_vel, ptl.g_vel)) - desc.vs;
            auto const alpha = std::acos(ptl.g_vel.x / (v + desc.vs));
            auto const g_ref = n * std::exp(-v * v / vth2) * std::pow(std::sin(alpha), zeta) / (M_PI * 2 * vth2 * std::sqrt(vth2) * Ab * Bz);
            REQUIRE(vdf.g0(ptl) == Approx{ g_ref }.epsilon(1e-10));
        }
    });

    if constexpr (dump_samples) {
        std::ofstream os{ "/Users/kyungguk/Downloads/RelativisticPartialShellVDF-homogeneous-anisotropic_shell.m" };
        os.setf(os.fixed);
        os.precision(20);
        static_assert(n_samples > 0);
        println(os, '{');
        for (unsigned long i = 0; i < particles.size() - 1; ++i) {
            println(os, "    ", particles[i], ", ");
        }
        println(os, "    ", particles.back());
        println(os, '}');
        os.close();
    }
}

TEST_CASE("Test LibPIC::RelativisticVDF::PartialShellVDF::Inhomogeneous::Maxwellian", "[LibPIC::RelativisticVDF::PartialShellVDF::Inhomogeneous::Maxwellian]")
{
    Real const O0 = 1, op = 4 * O0, c = op, beta = 0.1;
    Real const xi = .876, xiD1q1max = M_PI_2 * 0.8;
    long const q1min = -7, q1max = 15;
    auto const D1   = xiD1q1max / (xi * std::max(std::abs(q1min), std::abs(q1max)));
    auto const geo  = Geometry{ xi, D1, O0 };
    auto const desc = PartialShellPlasmaDesc({ { -O0, op }, 10, ShapeOrder::CIC }, beta);
    auto const vdf  = RelativisticPartialShellVDF(desc, geo, { q1min, q1max - q1min }, c);

    CHECK(serialize(desc) == serialize(vdf.plasma_desc()));

    // check equilibrium macro variables
    CHECK(vdf.Nrefcell_div_Ntotal() == Approx{ 1.0 / (q1max - q1min) }.epsilon(1e-10));

    for (long q1 = q1min; q1 <= q1max; ++q1) {
        CurviCoord const pos{ Real(q1) };
        auto const       eta = 1;

        auto const n0_ref = eta;
        auto const n0     = vdf.n0(pos);
        CHECK(*n0 == Approx{ n0_ref }.epsilon(1e-10));

        auto const nV0 = geo.cart_to_fac(vdf.nV0(pos), pos);
        CHECK(nV0.x == Approx{ 0 }.margin(1e-10));
        CHECK(nV0.y == Approx{ 0 }.margin(1e-10));
        CHECK(nV0.z == Approx{ 0 }.margin(1e-10));

        auto const nuv0 = geo.cart_to_fac(vdf.nuv0(pos), pos);
        auto const xs   = desc.vs / std::sqrt(beta);
        auto const xs2  = xs * xs;
        auto const Ab   = .5 * (xs * std::exp(-xs2) + 2 / M_2_SQRTPI * (.5 + xs2) * std::erfc(-xs));
        auto const T    = .5 * desc.beta / Ab * (xs * (2.5 + xs2) * std::exp(-xs2) + 2 / M_2_SQRTPI * (.75 + xs2 * (3 + xs2)) * std::erfc(-xs));
        auto const TT   = .5 * desc.beta * desc.beta / Ab * (xs * (33. / 4 + xs2 * (7 + xs2)) * std::exp(-xs2) + 2 / M_2_SQRTPI * (15. / 8 + xs2 * (45. / 4 + 15. / 2 * xs2 + xs2 * xs2)) * std::erfc(-xs));
        REQUIRE(nuv0.tt == Approx{ *n0 * (c * c + .5 * T) }.epsilon(1e-10));
        REQUIRE(nuv0.ts.x == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ts.y == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ts.z == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.xx == Approx{ *n0 / Real(3 + desc.zeta) * (T - .5 * TT / (c * c)) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.yy == Approx{ *n0 * (1 + .5 * desc.zeta) / Real(3 + desc.zeta) * (T - .5 * TT / (c * c)) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.zz == Approx{ *n0 * (1 + .5 * desc.zeta) / Real(3 + desc.zeta) * (T - .5 * TT / (c * c)) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.xy == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.yz == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.zx == Approx{ 0 }.margin(1e-10));
    }

    // sampling
    auto const n_samples = 100000U;
    auto const particles = vdf.emit(n_samples);

    static_assert(n_samples > 100);
    std::for_each_n(begin(particles), 100, [&](Particle const &ptl) {
        REQUIRE(ptl.psd.weight == 1);
        REQUIRE(ptl.psd.marker == vdf.f0(ptl));
        REQUIRE(ptl.psd.real_f == Approx{ ptl.psd.weight * ptl.psd.marker }.epsilon(1e-10));
        REQUIRE(vdf.real_f0(ptl) == vdf.f0(ptl));
        REQUIRE(ptl.gamma == Approx{ std::sqrt(1 + dot(ptl.g_vel, ptl.g_vel) / (c * c)) }.epsilon(1e-10));

        auto const n     = *vdf.n0(ptl.pos);
        auto const vth2  = beta;
        auto const v     = std::sqrt(dot(ptl.g_vel, ptl.g_vel));
        auto const f_ref = n * std::exp(-v * v / vth2) / (M_PI * 2 / M_2_SQRTPI * vth2 * std::sqrt(vth2));
        REQUIRE(vdf.f0(ptl) == Approx{ f_ref }.epsilon(1e-10));

        auto const marker_vth2 = vth2 * desc.marker_temp_ratio;
        auto const g_ref       = n * std::exp(-v * v / marker_vth2) / (M_PI * 2 / M_2_SQRTPI * marker_vth2 * std::sqrt(marker_vth2));
        REQUIRE(vdf.g0(ptl) == Approx{ g_ref }.epsilon(1e-10));
    });

    if constexpr (dump_samples) {
        std::ofstream os{ "/Users/kyungguk/Downloads/RelativisticPartialShellVDF-inhomogeneous-maxwellian.m" };
        os.setf(os.fixed);
        os.precision(20);
        static_assert(n_samples > 0);
        println(os, '{');
        for (unsigned long i = 0; i < particles.size() - 1; ++i) {
            println(os, "    ", particles[i], ", ");
        }
        println(os, "    ", particles.back());
        println(os, '}');
        os.close();
    }
}

TEST_CASE("Test LibPIC::RelativisticVDF::PartialShellVDF::Inhomogeneous::AnisotropicShell", "[LibPIC::RelativisticVDF::PartialShellVDF::Inhomogeneous::AnisotropicShell]")
{
    unsigned const zeta = 10;
    Real const     O0 = 1, op = 4 * O0, c = op, beta = 0.1, vs = 2;
    Real const     xi = .876, xiD1q1max = M_PI_2 * 0.8, psd_refresh_frequency = 0;
    long const     q1min = -7, q1max = 15;
    auto const     D1   = xiD1q1max / (xi * std::max(std::abs(q1min), std::abs(q1max)));
    auto const     geo  = Geometry{ xi, D1, O0 };
    auto const     desc = PartialShellPlasmaDesc({ { -O0, op }, 10, ShapeOrder::CIC, psd_refresh_frequency, delta_f, .001, 2.1 }, beta, zeta, vs);
    auto const     vdf  = RelativisticPartialShellVDF(desc, geo, { q1min, q1max - q1min }, c);

    auto const g_desc = PartialShellPlasmaDesc({ { -O0, op }, 10, ShapeOrder::CIC }, beta * desc.marker_temp_ratio, zeta, vs);
    auto const g_vdf  = RelativisticPartialShellVDF(g_desc, geo, { q1min, q1max - q1min }, c);

    CHECK(serialize(desc) == serialize(vdf.plasma_desc()));

    // check equilibrium macro variables
    CHECK(vdf.Nrefcell_div_Ntotal() == Approx{ 0.1108013193858871 }.epsilon(1e-10));

    std::array const etas{
        0.16070869470151475, 0.26703071776468223, 0.40485697458795117, 0.5642234318909161, 0.7267234017099303,
        0.868463747345441, 0.9654769492087792, 1., 0.9654769492087792, 0.868463747345441, 0.7267234017099303,
        0.5642234318909161, 0.40485697458795117, 0.26703071776468223, 0.16070869470151475, 0.08739009484781933,
        0.0423729187357854, 0.017993222556708512, 0.0065264207733627755, 0.001950941625094117, 0.00045560222393496486,
        0.00007636454203949487, 7.940057378694681e-6
    };
    long q1 = q1min;
    for (auto const &eta : etas) {
        CurviCoord const pos{ Real(q1++) };

        auto const n0_ref = eta;
        auto const n0     = vdf.n0(pos);
        CHECK(*n0 == Approx{ n0_ref }.epsilon(1e-10));

        auto const nV0 = geo.cart_to_fac(vdf.nV0(pos), pos);
        CHECK(nV0.x == Approx{ 0 }.margin(1e-10));
        CHECK(nV0.y == Approx{ 0 }.margin(1e-10));
        CHECK(nV0.z == Approx{ 0 }.margin(1e-10));

        auto const nuv0 = geo.cart_to_fac(vdf.nuv0(pos), pos);
        auto const xs   = desc.vs / std::sqrt(beta);
        auto const xs2  = xs * xs;
        auto const Ab   = .5 * (xs * std::exp(-xs2) + 2 / M_2_SQRTPI * (.5 + xs2) * std::erfc(-xs));
        auto const T    = .5 * desc.beta / Ab * (xs * (2.5 + xs2) * std::exp(-xs2) + 2 / M_2_SQRTPI * (.75 + xs2 * (3 + xs2)) * std::erfc(-xs));
        auto const TT   = .5 * desc.beta * desc.beta / Ab * (xs * (33. / 4 + xs2 * (7 + xs2)) * std::exp(-xs2) + 2 / M_2_SQRTPI * (15. / 8 + xs2 * (45. / 4 + 15. / 2 * xs2 + xs2 * xs2)) * std::erfc(-xs));
        REQUIRE(nuv0.tt == Approx{ *n0 * (c * c + .5 * T) }.epsilon(1e-10));
        REQUIRE(nuv0.ts.x == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ts.y == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ts.z == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.xx == Approx{ *n0 / Real(3 + desc.zeta) * (T - .5 * TT / (c * c)) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.yy == Approx{ *n0 * (1 + .5 * desc.zeta) / Real(3 + desc.zeta) * (T - .5 * TT / (c * c)) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.zz == Approx{ *n0 * (1 + .5 * desc.zeta) / Real(3 + desc.zeta) * (T - .5 * TT / (c * c)) }.epsilon(1e-10));
        REQUIRE(nuv0.ss.xy == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.yz == Approx{ 0 }.margin(1e-10));
        REQUIRE(nuv0.ss.zx == Approx{ 0 }.margin(1e-10));
    }

    // sampling
    auto const n_samples = 100000U;
    auto const particles = vdf.emit(n_samples);

    static_assert(n_samples > 100);
    std::for_each_n(begin(particles), 100, [&](Particle const &ptl) {
        REQUIRE(ptl.psd.weight == desc.initial_weight);
        REQUIRE(ptl.psd.marker == g_vdf.f0(ptl));
        REQUIRE(ptl.psd.real_f == Approx{ vdf.f0(ptl) + ptl.psd.weight * ptl.psd.marker }.epsilon(1e-10));
        REQUIRE(vdf.real_f0(ptl) == vdf.f0(ptl));
        REQUIRE(ptl.gamma == Approx{ std::sqrt(1 + dot(ptl.g_vel, ptl.g_vel) / (c * c)) }.epsilon(1e-10));

        auto const n = *vdf.n0(ptl.pos);
        {
            auto const vth2  = beta;
            auto const xs    = desc.vs / std::sqrt(vth2);
            auto const Ab    = .5 * (xs * std::exp(-xs * xs) + 2 / M_2_SQRTPI * (.5 + xs * xs) * std::erfc(-xs));
            auto const Bz    = 2 / M_2_SQRTPI * std::tgamma(1 + .5 * desc.zeta) / std::tgamma(1.5 + .5 * desc.zeta);
            auto const v     = std::sqrt(dot(ptl.g_vel, ptl.g_vel)) - desc.vs;
            auto const alpha = std::acos(ptl.g_vel.x / (v + desc.vs));
            auto const f_ref = n * std::exp(-v * v / vth2) * std::pow(std::sin(alpha), zeta) / (M_PI * 2 * vth2 * std::sqrt(vth2) * Ab * Bz);
            REQUIRE(vdf.f0(ptl) == Approx{ f_ref }.epsilon(1e-10));
        }
        {
            auto const vth2  = beta * desc.marker_temp_ratio;
            auto const xs    = desc.vs / std::sqrt(vth2);
            auto const Ab    = .5 * (xs * std::exp(-xs * xs) + 2 / M_2_SQRTPI * (.5 + xs * xs) * std::erfc(-xs));
            auto const Bz    = 2 / M_2_SQRTPI * std::tgamma(1 + .5 * desc.zeta) / std::tgamma(1.5 + .5 * desc.zeta);
            auto const v     = std::sqrt(dot(ptl.g_vel, ptl.g_vel)) - desc.vs;
            auto const alpha = std::acos(ptl.g_vel.x / (v + desc.vs));
            auto const g_ref = n * std::exp(-v * v / vth2) * std::pow(std::sin(alpha), zeta) / (M_PI * 2 * vth2 * std::sqrt(vth2) * Ab * Bz);
            REQUIRE(vdf.g0(ptl) == Approx{ g_ref }.epsilon(1e-10));
        }
    });

    if constexpr (dump_samples) {
        std::ofstream os{ "/Users/kyungguk/Downloads/RelativisticPartialShellVDF-inhomogeneous-anisotropic_shell.m" };
        os.setf(os.fixed);
        os.precision(20);
        static_assert(n_samples > 0);
        println(os, '{');
        for (unsigned long i = 0; i < particles.size() - 1; ++i) {
            println(os, "    ", particles[i], ", ");
        }
        println(os, "    ", particles.back());
        println(os, '}');
        os.close();
    }
}
#endif
