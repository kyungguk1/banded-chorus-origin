/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "catch2/catch.hpp"

#include <PIC/FourVector.h>
#include <PIC/Vector.h>

TEST_CASE("Test LibPIC::Vector", "[LibPIC::Vector]")
{
    {
        constexpr Vector v1{};
        constexpr bool   tf = v1.fold(true, [](bool lhs, auto rhs) {
            return lhs && rhs == 0;
          });
        CHECK(tf);

        constexpr Vector v2{ 1 };
        CHECK(v2.fold(true, [](bool lhs, auto rhs) {
            return lhs && rhs == 1;
        }));

        constexpr Vector v3{ 1, 2, 3 };
        CHECK((v3.x == 1 && v3.y == 2 && v3.z == 3));

        constexpr bool tf2 = std::addressof(v1) == std::addressof(+v1);
        CHECK(tf2);

        constexpr Vector v4 = -v3;
        CHECK((v4.x == -1 && v4.y == -2 && v4.z == -3));

        constexpr auto dot1 = dot(v1, v3);
        CHECK(dot1 == 0);
        constexpr auto dot2 = dot(v2, v3);
        CHECK(dot2 == 6);

        constexpr auto cross1 = cross(v1, v3);
        CHECK(cross1.fold(true, [](bool lhs, auto rhs) {
            return lhs && rhs == 0;
        }));
        constexpr auto cross2 = cross(Vector{ 4, 3, 8 }, v3);
        CHECK((cross2.x == -7 && cross2.y == -4 && cross2.z == 5));
    }

    {
        constexpr auto is_equal = [](Vector lhs, Vector rhs) {
            return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
        };

        constexpr Vector v1{ 1, 2, 3 };
        constexpr double x{ 1 };
        CHECK(is_equal(v1 + x, { 2, 3, 4 }));
        CHECK(is_equal(v1 - x, { 0, 1, 2 }));
        CHECK(is_equal(v1 * x, v1));
        CHECK(is_equal(v1 / x, v1));
        CHECK(is_equal(x + v1, v1 + x));
        CHECK(is_equal(x - v1, -(v1 - x)));
        CHECK(is_equal(x * v1, v1));
        CHECK(is_equal(x / v1, { x / 1, x / 2, x / 3 }));

        constexpr Vector v2 = v1 * 10.;
        CHECK(is_equal(v1 + v2, { 11, 22, 33 }));
        CHECK(is_equal(v2 - v1, { 9, 18, 27 }));
        CHECK(is_equal(v1 * v2, { 10, 40, 90 }));
        CHECK(is_equal(v2 / v1, { 10, 10, 10 }));
    }
}

TEST_CASE("Test LibPIC::FourVector", "[LibPIC::FourVector]")
{
    {
        constexpr FourVector v1{};
        CHECK((*v1.t == 0 && v1.s.x == 0 && v1.s.y == 0 && v1.s.z == 0));

        constexpr FourVector v2{ 1 };
        CHECK((*v2.t == 1 && v2.s.x == 1 && v2.s.y == 1 && v2.s.z == 1));

        constexpr FourVector v3{ 1, { 2, 3, 4 } };
        CHECK((*v3.t == 1 && v3.s.x == 2 && v3.s.y == 3 && v3.s.z == 4));

        constexpr bool tf2 = std::addressof(v1) == std::addressof(+v1);
        CHECK(tf2);

        constexpr FourVector v4 = -v3;
        CHECK((*v4.t == -1 && v4.s.x == -2 && v4.s.y == -3 && v4.s.z == -4));
    }

    {
        constexpr auto is_equal = [](FourVector const &lhs, FourVector const &rhs) {
            return *lhs.t == *rhs.t && lhs.s.x == rhs.s.x && lhs.s.y == rhs.s.y
                && lhs.s.z == rhs.s.z;
        };

        constexpr FourVector v1{ -1, { 1, 2, 3 } };
        constexpr double     x{ 1 };
        CHECK(is_equal(v1 + x, { 0, { 2, 3, 4 } }));
        CHECK(is_equal(v1 - x, { -2, { 0, 1, 2 } }));
        CHECK(is_equal(v1 * x, v1));
        CHECK(is_equal(v1 / x, v1));
        CHECK(is_equal(x + v1, v1 + x));
        CHECK(is_equal(x - v1, -(v1 - x)));
        CHECK(is_equal(x * v1, v1));
        CHECK(is_equal(x / v1, { x / -1, { x / 1, x / 2, x / 3 } }));

        constexpr FourVector v2 = v1 * 10.;
        CHECK(is_equal(v1 + v2, { -11, { 11, 22, 33 } }));
        CHECK(is_equal(v2 - v1, { -9, { 9, 18, 27 } }));
        CHECK(is_equal(v1 * v2, { 10, { 10, 40, 90 } }));
        CHECK(is_equal(v2 / v1, { 10, { 10, 10, 10 } }));
    }
}

namespace {
[[nodiscard]] constexpr bool operator==(std::complex<double> const &a, std::complex<double> const &b) noexcept
{
    return a.real() == b.real() && a.imag() == b.imag();
}
} // namespace
TEST_CASE("Test LibPIC::ComplexVector", "[LibPIC::ComplexVector]")
{
    {
        constexpr ComplexVector v1{};
        constexpr bool          tf
            = v1.fold(true, [](bool lhs, auto const &rhs) {
                  return lhs && rhs == 0;
              });
        CHECK(tf);

        constexpr ComplexVector v2{ 1i };
        CHECK(v2.fold(true, [](bool lhs, auto rhs) {
            return lhs && rhs == 1i;
        }));

        constexpr ComplexVector v3{ 1i, 2, 3i };
        CHECK((v3.x == 1i && v3.y == 2 && v3.z == 3i));

        constexpr bool tf2 = std::addressof(v1) == std::addressof(+v1);
        CHECK(tf2);

        ComplexVector const v4 = -v3;
        CHECK((v4.x == -1i && v4.y == -2 && v4.z == -3i));

        auto const dot1 = dot(v1, v3);
        CHECK(dot1 == 0.);
        auto const dot2 = dot(v2, v3);
        CHECK(dot2 == -4. + 2i);

        auto const cross1 = cross(v1, v3);
        CHECK(cross1.fold(true, [](bool lhs, auto rhs) {
            return lhs && rhs == 0;
        }));
        auto const cross2 = cross(ComplexVector{ 4, 3i, 8i }, v3);
        CHECK((cross2.x == -9. - 16.i && cross2.y == -8. - 12.i && cross2.z == 11. + 0.i));
    }

    {
        constexpr auto is_equal = [](ComplexVector lhs, ComplexVector rhs) {
            return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
        };

        constexpr ComplexVector             v1{ 1, 2i, 3 };
        constexpr ComplexVector::value_type z{ 1, 1 };
        CHECK(is_equal(v1 + z, { 2. + 1.i, 1. + 3.i, 4. + 1.i }));
        CHECK(is_equal(v1 - z, { 0. - 1.i, -1. + 1.i, 2. - 1.i }));
        CHECK(is_equal(v1 * z, { 1. + 1.i, -2. + 2.i, 3. + 3.i }));
        CHECK(is_equal(v1 / z, { 0.5 - 0.5i, 1. + 1.i, 1.5 - 1.5i }));
        CHECK(is_equal(z + v1, v1 + z));
        CHECK(is_equal(z - v1, -(v1 - z)));
        CHECK(is_equal(z * v1, v1 * z));
        CHECK(is_equal(z / v1, { z / v1.x, z / v1.y, z / v1.z }));

        ComplexVector const v2 = v1 * 10.;
        CHECK(is_equal(v1 + v2, { 11., 0. + 22.i, 33. }));
        CHECK(is_equal(v2 - v1, { 9., 0. + 18.i, 27. }));
        CHECK(is_equal(v1 * v2, { 10., -40. + 0.i, 90. }));
        CHECK(is_equal(v2 / v1, { 10., 10. + 0.i, 10. }));
    }
}
