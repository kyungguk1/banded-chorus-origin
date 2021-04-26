/*
 * Copyright (c) 2021, Kyungguk Min
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * - Redistributions of source code must retain the above copyright notice, this
 *  list of conditions and the following disclaimer.
 *
 * - Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
            throw std::domain_error{__PRETTY_FUNCTION__};

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
