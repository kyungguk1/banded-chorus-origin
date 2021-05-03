/*
 * Copyright (c) 2019, Kyungguk Min
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

#ifndef Current_h
#define Current_h

#include "../Geometry.h"
#include "../ParamSet.h"

HYBRID1D_BEGIN_NAMESPACE
class BField;
class EField;
class Lambda;
class Species;
class Gamma;

/// current density
///
class Current : public VectorGrid {
    VectorGrid tmp;

public:
    ParamSet const params;
    Geometry const geomtr;

public:
    virtual ~Current() = default;
    explicit Current(ParamSet const &);

    void reset() noexcept { this->fill(Vector{ 0 }); }
    void smooth() noexcept { _smooth(tmp, *this), this->swap(tmp); }

    virtual Current &operator+=(Species const &sp) noexcept;

    void advance(Lambda const &lambda, Gamma const &gamma, BField const &bfield,
                 EField const &efield, Real dt) noexcept;

private:
    static inline void _advance(Current &J, Lambda const &L, Gamma const &G, BField const &B,
                                EField const &E, Real dt) noexcept;
};

/// Γ
///
class Gamma : public Current {
    using Current::advance;
    using Current::smooth;

public:
    using Current::Current;
    Gamma &operator+=(Species const &sp) noexcept override;
};
HYBRID1D_END_NAMESPACE

#endif /* Current_h */
