/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC-config.h>

LIBPIC_BEGIN_NAMESPACE
/// Global real type
///
using Real = double;

/// Order of the shape function
///
enum ShapeOrder : long {
    CIC  = 1, //!< First order; cloud-in-cell scheme.
    TSC  = 2, //!< Second order; Triangular-shaped density cloud ssheme.
    _1st = 1, //!< Synonym for CIC.
    _2nd = 2, //!< Synonym for TSC.
    _3rd = 3, //!< 3rd order.
};

/// Flag indicating how to evolve the distribution function
///
enum ParticleScheme : bool {
    full_f  = 0, //!< Evolve full VDF.
    delta_f = 1, //!< Evolve the deviation of VDF from the initial.
};

/// Algorithm for electric field extrapolation.
///
enum Algorithm : long {
    PC,   //!< Using predictor-corrector by Kunz et al. (2014).
    CAMCL //!< Using CAM-CL by Matthew (1994).
};

/// Electron fluid closure
///
enum Closure : long {
    isothermal = 11, //!< gamma = 1/1.
    adiabatic  = 53  //!< gamma = 5/3.
};
LIBPIC_END_NAMESPACE
