//
//  Inputs.h
//  pic_1d
//
//  Created by KYUNGGUK MIN on 1/14/19.
//  Copyright © 2019 Kyungguk Min & Kaijun Liu. All rights reserved.
//

#ifndef Inputs_h
#define Inputs_h

/// simulation input parameters;
/// modify the variables under the `Input' namespace
/// consult "Predefined.h" header for symbol definitions and constants
///
namespace Input {
    //
    // MARK:- Housekeeping
    //

    /// number of worker threads for parallelization
    ///
    /// value `0' means serial update; value `n' means parallelization using n + 1 threads
    /// PtlDesc::Ncs*Nx must be divisible by n + 1
    ///
    constexpr unsigned number_of_worker_threads = 0;

    /// particle and interpolation order
    ///
    constexpr _ShapeOrder shape_order = CIC;

    /// number of source smoothings
    ///
    constexpr unsigned Nsmooths = 1;

    //
    // MARK: Global parameters
    //

    /// light speed
    ///
    constexpr Real c = 214.243;

    /// magnitude of uniform background magnetic field
    ///
    constexpr Real O0 = 1;

    /// angle in degrees between the x-axis and the uniform magnetic field direction.
    ///
    constexpr Real theta = 0;

    /// simulation grid size
    ///
    constexpr Real Dx = 0.2;

    /// number of grid points
    ///
    constexpr unsigned Nx = 1440;

    /// time step size
    ///
    constexpr Real dt = 0.01;

    /// number of time steps for inner loop
    /// total time step Nt = inner_Nt * outer_Nt
    /// simulation time t = dt*Nt
    ///
    constexpr unsigned inner_Nt = 20;

    /// number of time steps for outer loop
    /// total time step Nt = inner_Nt * outer_Nt
    /// simulation time t = dt*Nt
    ///
    constexpr unsigned outer_Nt = 5000;

    //
    // MARK: Particle Species Descriptions
    //
    namespace PtlDesc {
        /// number of particle species (or populations)
        ///
        constexpr unsigned Ns = 3;

        /// number of simulation particles per cell for individual populations
        ///
        constexpr std::array<unsigned, Ns> Ncs = {1000, 500, 500};

        /// species cyclotron frequencies for individual populations
        ///
        constexpr std::array<Real, Ns> Ocs = {1, 1, .25};

        /// species plasma frequencies for individual populations
        ///
        constexpr std::array<Real, Ns> ops = {47.9062, 207.716, 10.7122};

        /// parallel (w.r.t the background magnetic field direction)
        /// species betas for individual populations
        ///
        constexpr std::array<Real, Ns> betas = {0.15, 0.0094, 0.0001};

        /// species temperature anisotropies (T_perp/T_para) for individual populations
        ///
        constexpr std::array<Real, Ns> T2OT1s = {3, 1, 1};
    }

    //
    // MARK: Data Recording
    //

    /// a top-level directory to which outputs will be saved
    ///
    constexpr char working_directory[] = "./data";

    /// field and particle energy density recording frequency; in units of inner_Nt
    /// `0' means `not interested'
    ///
    constexpr unsigned energy_recording_frequency = 1;

    /// electric and magnetic field recording frequency
    ///
    constexpr unsigned field_recording_frequency = 2;

    /// kinetic ion moment recording frequency
    ///
    constexpr unsigned moment_recording_frequency = 10000;

    /// simulation particle recording frequency
    ///
    constexpr unsigned particle_recording_frequency = 10000;

    /// maximum number of particles to dump
    ///
    constexpr std::array<unsigned, PtlDesc::Ns> Ndumps = {1000, 900, 500};
}

#endif /* Inputs_h */
