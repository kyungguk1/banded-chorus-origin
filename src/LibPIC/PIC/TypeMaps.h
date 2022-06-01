/*
 * Copyright (c) 2021-2022, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC/CartCoord.h>
#include <PIC/Config.h>
#include <PIC/CurviCoord.h>
#include <PIC/FourTensor.h>
#include <PIC/FourVector.h>
#include <PIC/Particle.h>
#include <PIC/Predefined.h>
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
namespace pic_ver_1 = LIBPIC_NAMESPACE(1);

template <>
struct TypeMap<pic_ver_1::Scalar> {
    using type = pic_ver_1::Scalar;
    using root = pic_ver_1::Real;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == alignof(root),
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
template <>
struct TypeMap<pic_ver_1::Vector> {
    using type = pic_ver_1::Vector;
    using root = std::array<pic_ver_1::Real, 3>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == 8,
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>().realigned(alignof(type)); }
};
template <>
struct TypeMap<pic_ver_1::FourVector> {
    using type = pic_ver_1::FourVector;
    using root = std::array<pic_ver_1::Real, 4>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == 8,
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>().realigned(alignof(type)); }
};
template <>
struct TypeMap<pic_ver_1::Tensor> {
    using type = pic_ver_1::Tensor;
    using root = std::array<pic_ver_1::Real, 6>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == 8,
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>().realigned(alignof(type)); }
};
template <>
struct TypeMap<pic_ver_1::FourTensor> {
    using type = pic_ver_1::FourTensor;
    using root = std::array<pic_ver_1::Real, 10>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == 8,
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>().realigned(alignof(type)); }
};
template <>
struct TypeMap<pic_ver_1::CartCoord> {
    using type = pic_ver_1::CartCoord;
    using root = std::array<pic_ver_1::Real, 1>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == 8,
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
template <>
struct TypeMap<pic_ver_1::CurviCoord> {
    using type = pic_ver_1::CurviCoord;
    using root = std::array<pic_ver_1::Real, 1>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == 8,
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
template <>
struct TypeMap<pic_ver_1::Particle::PSD> {
    using type = pic_ver_1::Particle::PSD;
    using root = std::array<pic_ver_1::Real, 3>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == alignof(root),
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
template <>
struct TypeMap<pic_ver_1::Particle> {
    using type         = pic_ver_1::Particle;
    using equivalent_t = std::array<pic_ver_1::Real, 8>;
    static_assert(sizeof(type) == sizeof(equivalent_t) && alignof(type) == 8,
                  "Custom TypeMap: invalid type signature");
    static constexpr type v{};
    [[nodiscard]] auto    operator()() const
    {
        auto t = make_type(v.vel, v.pos, v.psd, v.id).realigned(alignof(type));
        if (t.extent().second != sizeof(type))
            throw std::domain_error{ __PRETTY_FUNCTION__ };
        return t;
    }
};
template <>
struct TypeMap<pic_ver_1::RelativisticParticle::PSD> {
    using type = pic_ver_1::RelativisticParticle::PSD;
    using root = std::array<pic_ver_1::Real, 3>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == alignof(root),
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
template <>
struct TypeMap<pic_ver_1::RelativisticParticle> {
    using type         = pic_ver_1::RelativisticParticle;
    using equivalent_t = std::array<pic_ver_1::Real, 9>;
    static_assert(sizeof(type) == sizeof(equivalent_t) && alignof(type) == 8,
                  "Custom TypeMap: invalid type signature");
    static constexpr type v{};
    [[nodiscard]] auto    operator()() const
    {
        auto t = make_type(v.g_vel, v.pos, v.psd, v.gamma, v.id).realigned(alignof(type));
        if (t.extent().second != sizeof(type))
            throw std::domain_error{ __PRETTY_FUNCTION__ };
        return t;
    }
};
} // namespace parallel

// hdf5 TypeMap Interfaces
//
namespace hdf5 {
namespace pic_ver_1 = LIBPIC_NAMESPACE(1);

template <>
struct TypeMap<pic_ver_1::Scalar> {
    using type = pic_ver_1::Scalar;
    using root = pic_ver_1::Real;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == alignof(root),
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
template <>
struct TypeMap<pic_ver_1::Vector> {
    using type = pic_ver_1::Vector;
    using root = std::array<pic_ver_1::Real, 3>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == 8,
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
template <>
struct TypeMap<pic_ver_1::FourVector> {
    using type = pic_ver_1::FourVector;
    using root = std::array<pic_ver_1::Real, 4>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == 8,
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
template <>
struct TypeMap<pic_ver_1::Tensor> {
    using type = pic_ver_1::Tensor;
    using root = std::array<pic_ver_1::Real, 6>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == 8,
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
template <>
struct TypeMap<pic_ver_1::FourTensor> {
    using type = pic_ver_1::FourTensor;
    using root = std::array<pic_ver_1::Real, 10>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == 8,
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
template <>
struct TypeMap<pic_ver_1::CartCoord> {
    using type = pic_ver_1::CartCoord;
    using root = std::array<pic_ver_1::Real, 1>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == 8,
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
template <>
struct TypeMap<pic_ver_1::CurviCoord> {
    using type = pic_ver_1::CurviCoord;
    using root = std::array<pic_ver_1::Real, 1>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == 8,
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
#if 1
template <>
struct TypeMap<pic_ver_1::Particle::PSD> {
    using type = pic_ver_1::Particle::PSD;
    using root = std::array<pic_ver_1::Real, 3>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == alignof(root),
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
template <>
struct TypeMap<pic_ver_1::Particle> {
    using type         = pic_ver_1::Particle;
    using equivalent_t = std::array<pic_ver_1::Real, 8>;
    static_assert(sizeof(type) == sizeof(equivalent_t) && alignof(type) == 8,
                  "Custom TypeMap: invalid type signature");
    static constexpr type v{};
    [[nodiscard]] auto    operator()() const
    {
        auto t = Type::compound(sizeof(type));
        t.insert("vel", HOFFSET(type, vel), make_type(v.vel));
        t.insert("pos", HOFFSET(type, pos), make_type(v.pos));
        t.insert("psd", HOFFSET(type, psd), make_type(v.psd));
        t.insert("id", HOFFSET(type, id), make_type(v.id));
        return t;
    }
};
template <>
struct TypeMap<pic_ver_1::RelativisticParticle::PSD> {
    using type = pic_ver_1::RelativisticParticle::PSD;
    using root = std::array<pic_ver_1::Real, 3>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == alignof(root),
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
template <>
struct TypeMap<pic_ver_1::RelativisticParticle> {
    using type         = pic_ver_1::RelativisticParticle;
    using equivalent_t = std::array<pic_ver_1::Real, 9>;
    static_assert(sizeof(type) == sizeof(equivalent_t) && alignof(type) == 8,
                  "Custom TypeMap: invalid type signature");
    static constexpr type v{};
    [[nodiscard]] auto    operator()() const
    {
        auto t = Type::compound(sizeof(type));
        t.insert("g_vel", HOFFSET(type, g_vel), make_type(v.g_vel));
        t.insert("pos", HOFFSET(type, pos), make_type(v.pos));
        t.insert("psd", HOFFSET(type, psd), make_type(v.psd));
        t.insert("gamma", HOFFSET(type, gamma), make_type(v.gamma));
        t.insert("id", HOFFSET(type, id), make_type(v.id));
        return t;
    }
};
#endif
} // namespace hdf5
