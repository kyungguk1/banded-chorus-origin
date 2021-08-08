/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#ifndef PIC_1D_TYPEMAPS_H
#define PIC_1D_TYPEMAPS_H

#include "./Particle.h"
#include "./Scalar.h"
#include "./Tensor.h"
#include "./Vector.h"

#include <HDF5Kit/HDF5Kit.h>
#include <ParallelKit/ParallelKit.h>
#include <array>
#include <stdexcept>
#include <type_traits>

// mpi TypeMap Interfaces
//
namespace parallel {
template <> struct TypeMap<P1D::Scalar> {
    using type = P1D::Scalar;
    using root = P1D::Real;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == alignof(root),
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
template <> struct TypeMap<P1D::Vector> {
    using type = P1D::Vector;
    using root = std::array<P1D::Real, 3>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == alignof(root),
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
template <> struct TypeMap<P1D::Tensor> {
    using type = P1D::Tensor;
    using root = std::array<P1D::Real, 6>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == alignof(root),
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
template <> struct TypeMap<P1D::Particle> {
    using type = P1D::Particle;
    static constexpr type v{};

    [[nodiscard]] auto operator()() const
    {
        auto t = make_type(v.vel, v.pos_x, v.f, v.w);
        if (t.extent().second != sizeof(type))
            throw std::domain_error{ __PRETTY_FUNCTION__ };

        return t;
    }
};
} // namespace parallel

// hdf5 TypeMap Interfaces
//
namespace hdf5 {
template <> struct TypeMap<P1D::Scalar> {
    using type = P1D::Scalar;
    using root = P1D::Real;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == alignof(root),
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
template <> struct TypeMap<P1D::Vector> {
    using type = P1D::Vector;
    using root = std::array<P1D::Real, 3>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == alignof(root),
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
template <> struct TypeMap<P1D::Tensor> {
    using type = P1D::Tensor;
    using root = std::array<P1D::Real, 6>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == alignof(root),
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
#if 0
template <> struct TypeMap<P1D::Particle> {
    using type         = P1D::Particle;
    using equivalent_t = std::array<P1D::Real, 6>;
    static_assert(sizeof(type) == sizeof(equivalent_t) && alignof(type) == alignof(equivalent_t));

    [[nodiscard]] auto operator()() const
    {
        auto t = Type::compound(sizeof(type));
        {
            static constexpr type v{};
            t.insert("vel", HOFFSET(type, vel), make_type(v.vel));
            t.insert("pos_x", HOFFSET(type, pos_x), make_type(v.pos_x));
            t.insert("f", HOFFSET(type, f), make_type(v.f));
            t.insert("w", HOFFSET(type, w), make_type(v.w));
        }
        return t;
    }
};
#endif
} // namespace hdf5

#endif // PIC_1D_TYPEMAPS_H
