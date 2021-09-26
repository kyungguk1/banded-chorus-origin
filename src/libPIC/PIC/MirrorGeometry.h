/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC/Config.h>
#include <PIC/MirrorBasis.h>
#include <PIC/MirrorCotrans.h>
#include <PIC/MirrorField.h>
#include <PIC/MirrorTransform.h>
#include <PIC/Predefined.h>
#include <PIC/Vector.h>

#include <cmath>

LIBPIC_BEGIN_NAMESPACE
/// Describes mirror field geometry
///
class MirrorGeometry
: public Detail::MirrorCotrans
, public Detail::MirrorField<MirrorGeometry>
, public Detail::MirrorBasis<MirrorGeometry>
, public Detail::MirrorTransform<MirrorGeometry> {
public:
    static constexpr Real inhomogeneity_xi_threshold = 1e-5;

    MirrorGeometry() noexcept;
    MirrorGeometry(Real xi, Vector const &D);

    /// Construct a MirrorGeometry object
    /// \param xi Mirror field inhomogeneity.
    /// \param D1 Grid scale factor along the parallel curvilinear coordinate.
    ///
    MirrorGeometry(Real xi, Real D1)
    : MirrorGeometry(xi, { D1, 1, 1 }) {}

    [[nodiscard]] Real xi() const noexcept { return m_xi; }
    // xi^2
    [[nodiscard]] Real xi2() const noexcept { return m_xi2; }
    [[nodiscard]] bool is_homogeneous() const noexcept { return m_homogeneous; }

    [[nodiscard]] Vector D() const noexcept { return m_D; }
    [[nodiscard]] Real   D1() const noexcept { return m_D.x; }
    [[nodiscard]] Real   D2() const noexcept { return m_D.y; }
    [[nodiscard]] Real   D3() const noexcept { return m_D.z; }

    // 1/D
    [[nodiscard]] Vector inv_D() const noexcept { return m_inv_D; }
    [[nodiscard]] Real   inv_D1() const noexcept { return m_inv_D.x; }
    [[nodiscard]] Real   inv_D2() const noexcept { return m_inv_D.y; }
    [[nodiscard]] Real   inv_D3() const noexcept { return m_inv_D.z; }

    // √g
    [[nodiscard]] Real sqrt_g() const noexcept { return m_sqrt_g; }
    // g = det(g_ij)
    [[nodiscard]] Real det_gij() const noexcept { return m_det_gij; }

    /// Check if the parallel curvilinear coordinate is within the valid range
    /// \param pos Curvilinear coordinate.
    ///
    [[nodiscard]] bool is_valid(CurviCoord const &pos) const noexcept { return std::abs(xi() * D1() * pos.q1) < M_PI_2; }

private:
    Vector m_D;
    Vector m_inv_D;
    Real   m_xi;
    Real   m_xi2;
    Real   m_sqrt_g;
    Real   m_det_gij;
    bool   m_homogeneous{};
};
LIBPIC_END_NAMESPACE
