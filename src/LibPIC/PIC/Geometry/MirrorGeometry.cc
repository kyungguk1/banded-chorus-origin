/*
 * Copyright (c) 2021-2022, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "MirrorGeometry.h"
#include <limits>
#include <stdexcept>

LIBPIC_NAMESPACE_BEGIN(1)
namespace {
constexpr auto quiet_nan = std::numeric_limits<Real>::quiet_NaN();
}
Detail::MirrorGeometry::MirrorGeometry() noexcept
: m_D{ quiet_nan }
, m_inv_D{ quiet_nan }
, m_xi{ quiet_nan }
, m_xi2{ quiet_nan }
, m_sqrt_g{ quiet_nan }
, m_det_gij{ quiet_nan }
{
}

Detail::MirrorGeometry::MirrorGeometry(Real const xi, Vector const &D)
: MirrorCotrans{ xi < inhomogeneity_xi_threshold }
, m_homogeneous{ xi < inhomogeneity_xi_threshold }
{
    if (xi < 0)
        throw std::invalid_argument{ std::string{ __PRETTY_FUNCTION__ } + " - negative xi" };
    if (D.x <= 0)
        throw std::invalid_argument{ std::string{ __PRETTY_FUNCTION__ } + " - non-positive D1" };
    if (D.y <= 0)
        throw std::invalid_argument{ std::string{ __PRETTY_FUNCTION__ } + " - non-positive D2" };
    if (D.z <= 0)
        throw std::invalid_argument{ std::string{ __PRETTY_FUNCTION__ } + " - non-positive D3" };

    m_D       = D;
    m_inv_D   = 1 / m_D;
    m_xi      = xi;
    m_xi2     = xi * xi;
    m_sqrt_g  = m_D.x * m_D.y * m_D.z;
    m_det_gij = m_sqrt_g * m_sqrt_g;
}
LIBPIC_NAMESPACE_END(1)
