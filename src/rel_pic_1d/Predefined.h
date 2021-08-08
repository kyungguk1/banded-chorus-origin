/*
 * Copyright (c) 2019, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#ifndef Predefined_h
#define Predefined_h

#include "./Macros.h"

PIC1D_BEGIN_NAMESPACE
/// Real type.
///
using Real = double;

/// Number of ghost cells.
///
constexpr long Pad = 3;

/// Order of the shape function.
///
enum ShapeOrder : long {
    CIC  = 1, //!<  First order; cloud-in-cell scheme.
    TSC  = 2, //!< Second order; Triangular-shaped density cloud sheme.
    _1st = 1, //!< Synonym of CIC.
    _2nd = 2, //!< Synonym of TSC.
    _3rd = 3  //!< 3rd order.
};

/// full-f versus delta-f.
///
enum ParticleScheme : bool { full_f = 0, delta_f = 1 };
PIC1D_END_NAMESPACE

#endif /* Predefined_h */
