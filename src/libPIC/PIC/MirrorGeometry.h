/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC/Config.h>
#include <PIC/MirrorCotrans.h>
#include <PIC/Predefined.h>
#include <PIC/Vector.h>

#include <cmath>

LIBPIC_BEGIN_NAMESPACE
class MirrorGeometry : public Detail::MirrorCotrans {
public:
    MirrorGeometry() noexcept = default;
    MirrorGeometry(Real xi, Real D1);

    [[nodiscard]] Real xi() const noexcept { return m_xi; }
    [[nodiscard]] Real xi2() const noexcept { return m_xi2; }
    [[nodiscard]] bool is_homogeneous() const noexcept { return m_homogeneous; }

    [[nodiscard]] Vector D() const noexcept { return m_D; }
    [[nodiscard]] Real   D1() const noexcept { return m_D.x; }
    [[nodiscard]] Real   D2() const noexcept { return m_D.y; }
    [[nodiscard]] Real   D3() const noexcept { return m_D.z; }

    [[nodiscard]] Real sqrt_g() const noexcept { return m_sqrt_g; }
    [[nodiscard]] Real det_gij() const noexcept { return m_det_gij; }

    [[nodiscard]] bool is_valid(CurviCoord const &pos) const noexcept { return std::abs(xi() * D1() * pos.q1) < M_PI_2; }

private:
    Vector m_D;
    Real   m_xi;
    Real   m_xi2;
    Real   m_sqrt_g;
    Real   m_det_gij;
    bool   m_homogeneous;
};
LIBPIC_END_NAMESPACE
