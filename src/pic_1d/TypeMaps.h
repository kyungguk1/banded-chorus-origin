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

#include "Utility/Particle.h"
#include "Utility/Scalar.h"
#include "Utility/Tensor.h"
#include "Utility/Vector.h"

#include <ParallelKit/ParallelKit.h>

// mpi TypeMap Interfaces
//
namespace parallel {
template <> struct TypeMap<P1D::Scalar> {
    using type = P1D::Scalar;
    using root = std::array<P1D::Real, 1>;
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
    using root = std::pair<std::array<P1D::Real, 4 /*{vx, vy, vz, x}*/>,
                           std::array<P1D::Real, 2 /*{f, w}*/>>;
    static_assert(sizeof(type) == sizeof(root) && alignof(type) == alignof(root),
                  "Custom TypeMap: invalid type signature");
    [[nodiscard]] auto operator()() const { return make_type<root>(); }
};
} // namespace parallel

#endif // PIC_1D_TYPEMAPS_H
