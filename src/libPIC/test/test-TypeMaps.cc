/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "catch2/catch.hpp"

#include <PIC/TypeMaps.h>
#include <memory>
#include <string>

TEST_CASE("Test libPIC::TypeMaps::ParallelKit", "[libPIC::TypeMaps::ParallelKit]")
{
    using parallel::make_type;

    try {
        using T      = Scalar;
        auto const t = make_type<T>();
        REQUIRE(!!t);
        CHECK(t.alignment() == alignof(T));
        CHECK(t.signature_size() == sizeof(T));

        auto [lb, extent] = t.extent();
        CHECK(lb == 0);
        CHECK(extent == sizeof(T));

        std::tie(lb, extent) = t.true_extent();
        CHECK(lb == 0);
        REQUIRE(extent == sizeof(T));
    } catch (std::exception const &e) {
        INFO(e.what())
        CHECK(false);
    }

    try {
        using T      = Vector;
        auto const t = make_type<T>();
        REQUIRE(!!t);
        CHECK(t.alignment() == alignof(T));
        CHECK(t.signature_size() == sizeof(T));

        auto [lb, extent] = t.extent();
        CHECK(lb == 0);
        CHECK(extent == sizeof(T));

        std::tie(lb, extent) = t.true_extent();
        CHECK(lb == 0);
        REQUIRE(extent == sizeof(T));
    } catch (std::exception const &e) {
        INFO(e.what())
        CHECK(false);
    }

    try {
        using T      = FourVector;
        auto const t = make_type<T>();
        REQUIRE(!!t);
        CHECK(t.alignment() == alignof(T));
        CHECK(t.signature_size() == sizeof(T));

        auto [lb, extent] = t.extent();
        CHECK(lb == 0);
        CHECK(extent == sizeof(T));

        std::tie(lb, extent) = t.true_extent();
        CHECK(lb == 0);
        REQUIRE(extent == sizeof(T));
    } catch (std::exception const &e) {
        INFO(e.what())
        CHECK(false);
    }

    try {
        using T      = Tensor;
        auto const t = make_type<T>();
        REQUIRE(!!t);
        CHECK(t.alignment() == alignof(T));
        CHECK(t.signature_size() == sizeof(T));

        auto [lb, extent] = t.extent();
        CHECK(lb == 0);
        CHECK(extent == sizeof(T));

        std::tie(lb, extent) = t.true_extent();
        CHECK(lb == 0);
        REQUIRE(extent == sizeof(T));
    } catch (std::exception const &e) {
        INFO(e.what())
        CHECK(false);
    }

    try {
        using T      = FourTensor;
        auto const t = make_type<T>();
        REQUIRE(!!t);
        CHECK(t.alignment() == alignof(T));
        CHECK(t.signature_size() == sizeof(T));

        auto [lb, extent] = t.extent();
        CHECK(lb == 0);
        CHECK(extent == sizeof(T));

        std::tie(lb, extent) = t.true_extent();
        CHECK(lb == 0);
        REQUIRE(extent == sizeof(T));
    } catch (std::exception const &e) {
        INFO(e.what())
        CHECK(false);
    }

    try {
        using T      = CartCoord;
        auto const t = make_type<T>();
        REQUIRE(!!t);
        CHECK(t.alignment() == alignof(T));
        CHECK(t.signature_size() == sizeof(T));

        auto [lb, extent] = t.extent();
        CHECK(lb == 0);
        CHECK(extent == sizeof(T));

        std::tie(lb, extent) = t.true_extent();
        CHECK(lb == 0);
        REQUIRE(extent == sizeof(T));
    } catch (std::exception const &e) {
        INFO(e.what())
        CHECK(false);
    }

    try {
        using T      = CurviCoord;
        auto const t = make_type<T>();
        REQUIRE(!!t);
        CHECK(t.alignment() == alignof(T));
        CHECK(t.signature_size() == sizeof(T));

        auto [lb, extent] = t.extent();
        CHECK(lb == 0);
        CHECK(extent == sizeof(T));

        std::tie(lb, extent) = t.true_extent();
        CHECK(lb == 0);
        REQUIRE(extent == sizeof(T));
    } catch (std::exception const &e) {
        INFO(e.what())
        CHECK(false);
    }

    try {
        using T      = Particle;
        auto const t = make_type<T>();
        REQUIRE(!!t);
        CHECK(t.alignment() == alignof(T));
        CHECK(t.signature_size() == sizeof(T));

        auto [lb, extent] = t.extent();
        CHECK(lb == 0);
        CHECK(extent == sizeof(T));

        std::tie(lb, extent) = t.true_extent();
        CHECK(lb == 0);
        REQUIRE(extent == sizeof(T));
    } catch (std::exception const &e) {
        INFO(e.what())
        CHECK(false);
    }

    try {
        using T      = RelativisticParticle;
        auto const t = make_type<T>();
        REQUIRE(!!t);
        CHECK(t.alignment() == alignof(T));
        CHECK(t.signature_size() == sizeof(T));

        auto [lb, extent] = t.extent();
        CHECK(lb == 0);
        CHECK(extent == sizeof(T));

        std::tie(lb, extent) = t.true_extent();
        CHECK(lb == 0);
        REQUIRE(extent == sizeof(T));
    } catch (std::exception const &e) {
        INFO(e.what())
        CHECK(false);
    }
}
TEST_CASE("Test libPIC::TypeMaps::HDF5Kit", "[libPIC::TypeMaps::HDF5Kit]")
{
    using hdf5::make_type;

    try {
        using T      = Scalar;
        auto const t = make_type<T>();
        REQUIRE(!!t);
        CHECK(t.size() == sizeof(T));
    } catch (std::exception const &e) {
        INFO("Exception thrown: " << e.what());
        REQUIRE(false);
    }

    try {
        using T      = Vector;
        auto const t = make_type<T>();
        REQUIRE(!!t);
        CHECK(t.size() == sizeof(T));
        CHECK(t.class_() == H5T_ARRAY);
        CHECK(H5Tequal(*t.super_(), H5T_NATIVE_DOUBLE));
    } catch (std::exception const &e) {
        INFO("Exception thrown: " << e.what());
        REQUIRE(false);
    }

    try {
        using T      = FourVector;
        auto const t = make_type<T>();
        REQUIRE(!!t);
        CHECK(t.size() == sizeof(T));
        CHECK(t.class_() == H5T_ARRAY);
        CHECK(H5Tequal(*t.super_(), H5T_NATIVE_DOUBLE));
    } catch (std::exception const &e) {
        INFO("Exception thrown: " << e.what());
        REQUIRE(false);
    }

    try {
        using T      = Tensor;
        auto const t = make_type<T>();
        REQUIRE(!!t);
        CHECK(t.size() == sizeof(T));
        CHECK(t.class_() == H5T_ARRAY);
        CHECK(H5Tequal(*t.super_(), H5T_NATIVE_DOUBLE));
    } catch (std::exception const &e) {
        INFO("Exception thrown: " << e.what());
        REQUIRE(false);
    }

    try {
        using T      = FourTensor;
        auto const t = make_type<T>();
        REQUIRE(!!t);
        CHECK(t.size() == sizeof(T));
        CHECK(t.class_() == H5T_ARRAY);
        CHECK(H5Tequal(*t.super_(), H5T_NATIVE_DOUBLE));
    } catch (std::exception const &e) {
        INFO("Exception thrown: " << e.what());
        REQUIRE(false);
    }

    try {
        using T      = CartCoord;
        auto const t = make_type<T>();
        REQUIRE(!!t);
        CHECK(t.size() == sizeof(T));
        CHECK(t.class_() == H5T_ARRAY);
        CHECK(H5Tequal(*t.super_(), H5T_NATIVE_DOUBLE));
    } catch (std::exception const &e) {
        INFO("Exception thrown: " << e.what());
        REQUIRE(false);
    }

    try {
        using T      = CurviCoord;
        auto const t = make_type<T>();
        REQUIRE(!!t);
        CHECK(t.size() == sizeof(T));
        CHECK(t.class_() == H5T_ARRAY);
        CHECK(H5Tequal(*t.super_(), H5T_NATIVE_DOUBLE));
    } catch (std::exception const &e) {
        INFO("Exception thrown: " << e.what());
        REQUIRE(false);
    }

    try {
        using T      = Particle;
        auto const t = make_type<T>();
        REQUIRE(!!t);
        CHECK(t.size() == sizeof(T));
        CHECK(t.class_() == H5T_COMPOUND);

        auto const n_members = H5Tget_nmembers(*t);
        REQUIRE(n_members == 4);

        using char_ptr = std::unique_ptr<char, void (*)(void *)>;
        char_ptr name{ nullptr, &free };

        CHECK(H5Tget_member_class(*t, 0) == H5T_ARRAY); // vel
        name = char_ptr{ H5Tget_member_name(*t, 0), &free };
        CHECK((!!name && std::string{ "vel" } == name.get()));

        CHECK(H5Tget_member_class(*t, 1) == H5T_ARRAY); // pos
        name = char_ptr{ H5Tget_member_name(*t, 1), &free };
        CHECK((!!name && std::string{ "pos" } == name.get()));

        CHECK(H5Tget_member_class(*t, 2) == H5T_ARRAY); // psd
        name = char_ptr{ H5Tget_member_name(*t, 2), &free };
        CHECK((!!name && std::string{ "psd" } == name.get()));

        CHECK(H5Tget_member_class(*t, 3) == H5T_INTEGER); // id
        name = char_ptr{ H5Tget_member_name(*t, 3), &free };
        CHECK((!!name && std::string{ "id" } == name.get()));
    } catch (std::exception const &e) {
        INFO("Exception thrown: " << e.what());
        REQUIRE(false);
    }

    try {
        using T      = RelativisticParticle;
        auto const t = make_type<T>();
        REQUIRE(!!t);
        CHECK(t.size() == sizeof(T));
        CHECK(t.class_() == H5T_COMPOUND);

        auto const n_members = H5Tget_nmembers(*t);
        REQUIRE(n_members == 5);

        using char_ptr = std::unique_ptr<char, void (*)(void *)>;
        char_ptr name{ nullptr, &free };

        CHECK(H5Tget_member_class(*t, 0) == H5T_ARRAY); // vel
        name = char_ptr{ H5Tget_member_name(*t, 0), &free };
        CHECK((!!name && std::string{ "g_vel" } == name.get()));

        CHECK(H5Tget_member_class(*t, 1) == H5T_ARRAY); // pos
        name = char_ptr{ H5Tget_member_name(*t, 1), &free };
        CHECK((!!name && std::string{ "pos" } == name.get()));

        CHECK(H5Tget_member_class(*t, 2) == H5T_ARRAY); // psd
        name = char_ptr{ H5Tget_member_name(*t, 2), &free };
        CHECK((!!name && std::string{ "psd" } == name.get()));

        CHECK(H5Tget_member_class(*t, 3) == H5T_FLOAT); // gamma
        name = char_ptr{ H5Tget_member_name(*t, 3), &free };
        CHECK((!!name && std::string{ "gamma" } == name.get()));

        CHECK(H5Tget_member_class(*t, 4) == H5T_INTEGER); // id
        name = char_ptr{ H5Tget_member_name(*t, 4), &free };
        CHECK((!!name && std::string{ "id" } == name.get()));
    } catch (std::exception const &e) {
        INFO("Exception thrown: " << e.what());
        REQUIRE(false);
    }
}
