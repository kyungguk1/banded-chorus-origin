/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "catch2/catch.hpp"

#include <PIC/PlasmaDesc.h>
#include <exception>

TEST_CASE("Test libPIC::PlasmaDesc", "[libPIC::PlasmaDesc]")
{
    constexpr auto desc1 = PlasmaDesc{ 1, 2, 3 };
    CHECK(desc1.Oc == 1);
    CHECK(desc1.op == 2);
    CHECK(desc1.number_of_source_smoothings == 3);

    constexpr auto desc2 = PlasmaDesc{ 1, 2 };
    CHECK(desc2.Oc == 1);
    CHECK(desc2.op == 2);
    CHECK(desc2.number_of_source_smoothings == 0);

    constexpr auto s1 = serialize(desc1);
    constexpr auto s2 = serialize(desc2);
    CHECK(s1 == s2);
    CHECK(desc1 == desc2);
    CHECK(std::get<0>(s1) == desc1.Oc);
    CHECK(std::get<1>(s1) == desc1.op);

    CHECK_THROWS_AS(PlasmaDesc(0, 1), std::exception);
    CHECK_THROWS_AS(PlasmaDesc(-1, 0), std::exception);
    CHECK_THROWS_AS(PlasmaDesc(-1, -1), std::exception);
}

TEST_CASE("Test libPIC::eFluidDesc", "[libPIC::eFluidDesc]")
{
    constexpr auto base1 = PlasmaDesc{ 1, 2, 3 };
    constexpr auto desc1 = eFluidDesc(base1, 1.1, Closure::isothermal);
    CHECK(desc1 == base1);
    CHECK(desc1.beta == 1.1);
    CHECK(desc1.gamma == 1);

    constexpr auto s1 = serialize(desc1);
    CHECK(std::get<2>(s1) == desc1.beta);
    CHECK(std::get<3>(s1) == desc1.gamma);

    constexpr auto base2 = PlasmaDesc{ 1, 2 };
    constexpr auto desc2 = eFluidDesc(base2);
    CHECK(desc2 == base2);
    CHECK(desc2.beta == 0);
    CHECK(desc2.gamma == 5. / 3);

    CHECK_THROWS_AS(eFluidDesc(base1, -1), std::exception);
}

TEST_CASE("Test libPIC::ColdPlasmaDesc", "[libPIC::ColdPlasmaDesc]")
{
    constexpr auto base1 = PlasmaDesc{ 1, 2, 3 };
    constexpr auto desc1 = ColdPlasmaDesc(base1);
    CHECK(desc1 == base1);

    constexpr auto base2 = PlasmaDesc{ 1, 2 };
    constexpr auto desc2 = ColdPlasmaDesc(base2);
    CHECK(desc2 == base2);
}

TEST_CASE("Test libPIC::KineticPlasmaDesc", "[libPIC::KineticPlasmaDesc]")
{
    constexpr auto base1 = PlasmaDesc{ 1, 2, 3 };
    constexpr auto desc1 = KineticPlasmaDesc(base1, 100, ShapeOrder::CIC);
    CHECK(desc1 == base1);
    CHECK(desc1.Nc == 100);
    CHECK(desc1.shape_order == 1);
    CHECK(desc1.scheme == ParticleScheme::full_f);
    CHECK(desc1.initial_weight == 0);
    CHECK(desc1.marker_temp_ratio == 1);

    constexpr auto base2 = PlasmaDesc{ 1, 2 };
    constexpr auto desc2
        = KineticPlasmaDesc(base2, 10, ShapeOrder::_3rd, ParticleScheme::delta_f, .1, 2);
    CHECK(desc2 == base2);
    CHECK(desc2.shape_order == ShapeOrder::_3rd);
    CHECK(desc2.scheme == ParticleScheme::delta_f);
    CHECK(desc2.initial_weight == .1);
    CHECK(desc2.marker_temp_ratio == 2);

    constexpr auto s1 = serialize(desc1);
    CHECK(std::get<2>(s1) == desc1.Nc);
    CHECK(std::get<3>(s1) == desc1.scheme);
    CHECK(desc1 == desc1);

    CHECK_THROWS_AS(KineticPlasmaDesc(base1, 0, ShapeOrder::TSC), std::exception);
    CHECK_THROWS_AS(KineticPlasmaDesc(base1, 0, ShapeOrder::TSC, ParticleScheme::delta_f, -1),
                    std::exception);
}

bool operator==(Vector const &a, Vector const &b)
{
    return a.x == b.x && a.y == b.y && a.z == b.z;
}
bool operator==(CurviCoord const &a, CurviCoord const &b)
{
    return a.q1 == b.q1;
}
TEST_CASE("Test libPIC::TestParticleDesc", "[libPIC::TestParticleDesc]")
{
    constexpr unsigned                      Nptls = 2;
    constexpr auto                          base1 = KineticPlasmaDesc{ { 1, 2, 3 }, 100, ShapeOrder::CIC };
    constexpr std::array<Vector, Nptls>     vel   = { Vector{ 1, 2, 3 }, { 4, 5, 6 } };
    constexpr std::array<CurviCoord, Nptls> pos   = { CurviCoord{ 1 }, CurviCoord{ 2 } };
    constexpr auto                          desc1 = TestParticleDesc<Nptls>(base1, vel, pos);
    CHECK_FALSE(desc1 == base1);
    CHECK(desc1.op == 0);
    CHECK(desc1.Nc == 0);
    CHECK(desc1.number_of_source_smoothings == 0);
    CHECK(desc1.number_of_test_particles == Nptls);
    CHECK(desc1.vel == vel);
    CHECK(desc1.pos == pos);

    (void)serialize(desc1);
    CHECK(desc1 == desc1);
}

TEST_CASE("Test libPIC::BiMaxPlasmaDesc", "[libPIC::BiMaxPlasmaDesc]")
{
    constexpr auto base1 = KineticPlasmaDesc{ { 1, 2, 3 }, 100, ShapeOrder::CIC };
    constexpr auto desc1 = BiMaxPlasmaDesc(base1, 1, 2);
    CHECK(desc1 == base1);
    CHECK(desc1.beta1 == 1);
    CHECK(desc1.T2_T1 == 2);

    constexpr auto base2 = KineticPlasmaDesc{ { 1, 2 }, 10, ShapeOrder::_3rd, ParticleScheme::delta_f };
    constexpr auto desc2 = BiMaxPlasmaDesc(base2, .1);
    CHECK(desc2 == base2);
    CHECK(desc2.beta1 == .1);
    CHECK(desc2.T2_T1 == 1);

    constexpr auto s1 = serialize(desc1);
    CHECK(std::get<6>(s1) == desc1.beta1);
    CHECK(std::get<7>(s1) == desc1.T2_T1);
    CHECK(desc1 == desc1);

    CHECK_THROWS_AS(BiMaxPlasmaDesc(base1, 0), std::exception);
    CHECK_THROWS_AS(BiMaxPlasmaDesc(base1, -1), std::exception);
}

TEST_CASE("Test libPIC::LossconePlasmaDesc", "[libPIC::LossconePlasmaDesc]")
{
    constexpr auto base1 = BiMaxPlasmaDesc({ { 1, 2, 3 }, 100, ShapeOrder::CIC }, 1, 2);
    constexpr auto desc1 = LossconePlasmaDesc(base1, .3);
    CHECK(desc1 == base1);
    CHECK(desc1.beta == .3);
    constexpr auto desc2 = LossconePlasmaDesc(base1);
    CHECK(desc2 == base1);
    CHECK(desc2.beta == 0);
    CHECK_THROWS_AS(LossconePlasmaDesc(base1, -1), std::exception);

    constexpr auto base2 = static_cast<KineticPlasmaDesc const &>(base1);
    constexpr auto desc3 = LossconePlasmaDesc(base2, base1.beta1);
    CHECK(desc3 == base2);
    CHECK(desc3.beta1 == 1);
    CHECK(desc3.T2_T1 == 1);
    CHECK(desc3.beta == 0);
    constexpr auto desc4 = LossconePlasmaDesc(base2, base1.beta1, 2);
    CHECK(desc4 == base2);
    CHECK(desc4.beta1 == 1);
    CHECK(desc4.T2_T1 == 2);
    CHECK(desc4.beta == 0);
    constexpr auto desc5 = LossconePlasmaDesc(base2, base1.beta1, 2, .1);
    CHECK(desc5 == base2);
    CHECK(desc5.beta1 == 1);
    CHECK(desc5.T2_T1 == 2 * (1 + .1));
    CHECK(desc5.beta == .1);
    constexpr auto desc6 = LossconePlasmaDesc(base2, base1.beta1, base1.T2_T1, 0);
    CHECK(desc6 == base1);

    CHECK_THROWS_AS(LossconePlasmaDesc(base2, .1, -1), std::exception);
    CHECK_THROWS_AS(LossconePlasmaDesc(base2, .1, 1, -1), std::exception);
}
