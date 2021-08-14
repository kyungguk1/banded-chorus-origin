/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <Utility/Scalar.h>
#include <Utility/Tensor.h>
#include <Utility/Vector.h>
#include <VDF/Particle.h>

#include <HDF5Kit/HDF5Kit.h>
#include <ParallelKit/ParallelKit.h>
#include <array>
#include <stdexcept>
#include <type_traits>

// mpi TypeMap Interfaces
//
namespace parallel {
template <> struct TypeMap<common::Scalar> {
    using type = common::Scalar;
    using root = common::Real;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == alignof(root),
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
template <> struct TypeMap<common::Vector> {
    using type = common::Vector;
    using root = std::array<common::Real, 4>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == 32,
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>().realigned(alignof(type)); }
};
template <> struct TypeMap<common::Tensor> {
    using type = common::Tensor;
    using root = std::array<common::Real, 8>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == 32,
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>().realigned(alignof(type)); }
};
template <> struct TypeMap<common::Particle::PSD> {
    using type = common::Particle::PSD;
    using root = std::array<common::Real, 2>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == alignof(root),
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
template <> struct TypeMap<common::Particle> {
    using type         = common::Particle;
    using equivalent_t = std::array<common::Real, 8>;
    static_assert(sizeof(type) == sizeof(equivalent_t) && alignof(type) == 32,
                  "Custom TypeMap: invalid type signature");
    static constexpr type v{};
    [[nodiscard]] auto    operator()() const
    {
        auto t = make_type(v.vel, v.pos_x, v.psd, v.id).realigned(alignof(type));
        if (t.extent().second != sizeof(type))
            throw std::domain_error{ __PRETTY_FUNCTION__ };
        return t;
    }
};
template <> struct TypeMap<common::RelativisticParticle> {
    using type         = common::RelativisticParticle;
    using equivalent_t = std::array<common::Real, 8>;
    static_assert(sizeof(type) == sizeof(equivalent_t) && alignof(type) == 32,
                  "Custom TypeMap: invalid type signature");
    static constexpr type v{};
    [[nodiscard]] auto    operator()() const
    {
        auto t = make_type(v.g_vel, v.pos_x, v.psd, v.gamma).realigned(alignof(type));
        if (t.extent().second != sizeof(type))
            throw std::domain_error{ __PRETTY_FUNCTION__ };
        return t;
    }
};
} // namespace parallel

// hdf5 TypeMap Interfaces
//
namespace hdf5 {
template <> struct TypeMap<common::Scalar> {
    using type = common::Scalar;
    using root = common::Real;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == alignof(root),
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
template <> struct TypeMap<common::Vector> {
    using type = common::Vector;
    using root = std::array<common::Real, 4>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == 32,
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
template <> struct TypeMap<common::Tensor> {
    using type = common::Tensor;
    using root = std::array<common::Real, 8>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == 32,
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
#if 1
template <> struct TypeMap<common::Particle::PSD> {
    using type = common::Particle::PSD;
    using root = std::array<common::Real, 2>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == alignof(root),
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
template <> struct TypeMap<common::Particle> {
    using type         = common::Particle;
    using equivalent_t = std::array<common::Real, 8>;
    static_assert(sizeof(type) == sizeof(equivalent_t) && alignof(type) == 32,
                  "Custom TypeMap: invalid type signature");
    static constexpr type v{};
    [[nodiscard]] auto    operator()() const
    {
        auto t = Type::compound(sizeof(type));
        t.insert("vel", HOFFSET(type, vel), make_type(v.vel));
        t.insert("pos_x", HOFFSET(type, pos_x), make_type(v.pos_x));
        t.insert("psd", HOFFSET(type, psd), make_type(v.psd));
        t.insert("id", HOFFSET(type, id), make_type(v.id));
        return t;
    }
};
template <> struct TypeMap<common::RelativisticParticle> {
    using type         = common::RelativisticParticle;
    using equivalent_t = std::array<common::Real, 8>;
    static_assert(sizeof(type) == sizeof(equivalent_t) && alignof(type) == 32,
                  "Custom TypeMap: invalid type signature");
    static constexpr type v{};
    [[nodiscard]] auto    operator()() const
    {
        auto t = Type::compound(sizeof(type));
        t.insert("g_vel", HOFFSET(type, g_vel), make_type(v.g_vel));
        t.insert("pos_x", HOFFSET(type, pos_x), make_type(v.pos_x));
        t.insert("psd", HOFFSET(type, psd), make_type(v.psd));
        t.insert("gamma", HOFFSET(type, gamma), make_type(v.gamma));
        return t;
    }
};
#endif
} // namespace hdf5