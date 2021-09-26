/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "MirrorGeometry.h"
#include <stdexcept>

LIBPIC_BEGIN_NAMESPACE
MirrorGeometry::MirrorGeometry(Real const xi, Vector const &D)
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
LIBPIC_END_NAMESPACE
