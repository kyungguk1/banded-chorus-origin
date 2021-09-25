/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "MirrorGeometry.h"
#include <stdexcept>

LIBPIC_BEGIN_NAMESPACE
MirrorGeometry::MirrorGeometry(Real const xi, Real const D1)
: MirrorCotrans{ xi < 1e-5 }
{
    if (xi < 0)
        throw std::invalid_argument{ std::string{ __PRETTY_FUNCTION__ } + " - negative xi" };
    if (D1 <= 0)
        throw std::invalid_argument{ std::string{ __PRETTY_FUNCTION__ } + " - non-positive D1" };

    m_D           = { D1, 1, 1 };
    m_xi          = xi;
    m_xi2         = xi * xi;
    m_sqrt_g      = m_D.x * m_D.y * m_D.z;
    m_det_gij     = m_sqrt_g * m_sqrt_g;
    m_homogeneous = m_xi2 < 1e-10;
}
LIBPIC_END_NAMESPACE
