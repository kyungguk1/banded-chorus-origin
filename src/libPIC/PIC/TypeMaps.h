/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC/FourTensor.h>
#include <PIC/FourVector.h>
#include <PIC/Particle.h>
#include <PIC/RelativisticParticle.h>
#include <PIC/Scalar.h>
#include <PIC/Tensor.h>
#include <PIC/Vector.h>

#include <HDF5Kit/HDF5Kit.h>
#include <ParallelKit/ParallelKit.h>
#include <array>
#include <stdexcept>
#include <type_traits>

// mpi TypeMap Interfaces
//
namespace parallel {
template <>
struct TypeMap<PIC::Scalar> {
    using type = PIC::Scalar;
    using root = PIC::Real;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == alignof(root),
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
template <>
struct TypeMap<PIC::Vector> {
    using type = PIC::Vector;
    using root = std::array<PIC::Real, 3>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == 8,
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>().realigned(alignof(type)); }
};
template <>
struct TypeMap<PIC::FourVector> {
    using type = PIC::FourVector;
    using root = std::array<PIC::Real, 4>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == 8,
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>().realigned(alignof(type)); }
};
template <>
struct TypeMap<PIC::Tensor> {
    using type = PIC::Tensor;
    using root = std::array<PIC::Real, 6>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == 8,
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>().realigned(alignof(type)); }
};
template <>
struct TypeMap<PIC::FourTensor> {
    using type = PIC::FourTensor;
    using root = std::array<PIC::Real, 10>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == 8,
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>().realigned(alignof(type)); }
};
template <>
struct TypeMap<PIC::Particle::PSD> {
    using type = PIC::Particle::PSD;
    using root = std::array<PIC::Real, 3>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == alignof(root),
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
template <>
struct TypeMap<PIC::Particle> {
    using type         = PIC::Particle;
    using equivalent_t = std::array<PIC::Real, 8>;
    static_assert(sizeof(type) == sizeof(equivalent_t) && alignof(type) == 8,
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
template <>
struct TypeMap<PIC::RelativisticParticle::PSD> {
    using type = PIC::RelativisticParticle::PSD;
    using root = std::array<PIC::Real, 3>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == alignof(root),
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
template <>
struct TypeMap<PIC::RelativisticParticle> {
    using type         = PIC::RelativisticParticle;
    using equivalent_t = std::array<PIC::Real, 9>;
    static_assert(sizeof(type) == sizeof(equivalent_t) && alignof(type) == 8,
                  "Custom TypeMap: invalid type signature");
    static constexpr type v{};
    [[nodiscard]] auto    operator()() const
    {
        auto t = make_type(v.g_vel, v.pos_x, v.psd, v.gamma, v.id).realigned(alignof(type));
        if (t.extent().second != sizeof(type))
            throw std::domain_error{ __PRETTY_FUNCTION__ };
        return t;
    }
};
} // namespace parallel

// hdf5 TypeMap Interfaces
//
namespace hdf5 {
template <>
struct TypeMap<PIC::Scalar> {
    using type = PIC::Scalar;
    using root = PIC::Real;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == alignof(root),
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
template <>
struct TypeMap<PIC::Vector> {
    using type = PIC::Vector;
    using root = std::array<PIC::Real, 3>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == 8,
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
template <>
struct TypeMap<PIC::FourVector> {
    using type = PIC::FourVector;
    using root = std::array<PIC::Real, 4>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == 8,
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
template <>
struct TypeMap<PIC::Tensor> {
    using type = PIC::Tensor;
    using root = std::array<PIC::Real, 6>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == 8,
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
template <>
struct TypeMap<PIC::FourTensor> {
    using type = PIC::FourTensor;
    using root = std::array<PIC::Real, 10>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == 8,
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
#if 1
template <>
struct TypeMap<PIC::Particle::PSD> {
    using type = PIC::Particle::PSD;
    using root = std::array<PIC::Real, 3>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == alignof(root),
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
template <>
struct TypeMap<PIC::Particle> {
    using type         = PIC::Particle;
    using equivalent_t = std::array<PIC::Real, 8>;
    static_assert(sizeof(type) == sizeof(equivalent_t) && alignof(type) == 8,
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
template <>
struct TypeMap<PIC::RelativisticParticle::PSD> {
    using type = PIC::RelativisticParticle::PSD;
    using root = std::array<PIC::Real, 3>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == alignof(root),
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
template <>
struct TypeMap<PIC::RelativisticParticle> {
    using type         = PIC::RelativisticParticle;
    using equivalent_t = std::array<PIC::Real, 9>;
    static_assert(sizeof(type) == sizeof(equivalent_t) && alignof(type) == 8,
                  "Custom TypeMap: invalid type signature");
    static constexpr type v{};
    [[nodiscard]] auto    operator()() const
    {
        auto t = Type::compound(sizeof(type));
        t.insert("g_vel", HOFFSET(type, g_vel), make_type(v.g_vel));
        t.insert("pos_x", HOFFSET(type, pos_x), make_type(v.pos_x));
        t.insert("psd", HOFFSET(type, psd), make_type(v.psd));
        t.insert("gamma", HOFFSET(type, gamma), make_type(v.gamma));
        t.insert("id", HOFFSET(type, id), make_type(v.id));
        return t;
    }
};
#endif
} // namespace hdf5
