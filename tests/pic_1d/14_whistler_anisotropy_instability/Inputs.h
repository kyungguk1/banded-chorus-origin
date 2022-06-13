/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

/// Simulation input parameters
/// \details
/// Modify the variables under the `Input' struct.
///
/// Consult/search headers under libPIC for unknown types/symbols.
///
struct Input {
    //
    // MARK:- Environment
    //

    /// number of subdomains for domain decomposition (positive integer)
    ///
    /// Nx must be divisible by this number
    ///
    static constexpr unsigned number_of_subdomains = 2;

    /// number of subdomain clones on which evenly divided particles are assigned and updated (positive integer)
    ///
    static constexpr unsigned number_of_distributed_particle_subdomain_clones = 3;

    /// number of worker threads to spawn for parallelization
    ///
    /// value `0' means serial update; value `n' means parallelization using n + 1 threads
    /// n + 1 must be divisible by number_of_subdomains * number_of_distributed_particle_subdomain_clones
    ///
    static constexpr unsigned number_of_worker_threads = 2 * number_of_subdomains * number_of_distributed_particle_subdomain_clones - 1;

    /// flag to suppress longitudinal and/or transverse components of the field fluctuations
    ///
    /// setting both to true will zero out all field fluctuations
    ///
    static constexpr bool should_neglect_transverse_component   = false;
    static constexpr bool should_neglect_longitudinal_component = false;

    /// particle boundary condition
    ///
    static constexpr BC particle_boundary_condition = BC::reflecting;

    /// if set, randomize the gyro-phase of reflected particles
    ///
    static constexpr bool should_randomize_gyrophase_of_reflecting_particles = false;

    /// wave masking function
    ///
    /// the first argument is masking inset, i.e., the number of grid points through which waves are gradually damped
    /// the second argument is the masking coefficients, zero being no masking at all and one being 0 to 100% masking within the masking inset
    ///
    /// see `docs/boundary_condition.nb` for how the phase retardation and amplitude damping work
    ///
    static constexpr MaskingFunction phase_retardation{ 50, 0.5 };
    static constexpr MaskingFunction amplitude_damping{ 50, 0.2 };

    //
    // MARK: Global parameters
    //

    /// light speed
    ///
    static constexpr Real c = 4;

    /// magnitude of equatorial background magnetic field
    ///
    static constexpr Real O0 = 1;

    /// inhomogeneity parameter, ξ
    /// the field variation at the central field line is given by B/B_eq = 1 + (ξ*x)^2,
    /// where x is the coordinate along the axis of mirror field symmetry
    ///
    static constexpr Real xi = 0.00260417;

    /// simulation grid size at the equator
    ///
    static constexpr Real Dx = 0.8;

    /// number of grid points
    ///
    static constexpr unsigned Nx = 960;

    /// time step size
    ///
    static constexpr Real dt = 0.02;

    /// number of time steps for inner loop
    /// total time step Nt = inner_Nt * outer_Nt
    /// simulation time t = dt*Nt
    ///
    static constexpr unsigned inner_Nt = 5;

    /// number of time steps for outer loop
    /// total time step Nt = inner_Nt * outer_Nt
    /// simulation time t = dt*Nt
    ///
    /// This is configurable through the command line option, `--outer_Nt=[0-9].*`
    ///
    static constexpr unsigned outer_Nt = 4000;

    //
    // MARK: Plasma Species Descriptions
    //

    /// kinetic plasma descriptors
    ///
    static constexpr auto part_descs
        = std::make_tuple(BiMaxPlasmaDesc({ { -O0, c / M_SQRT2, 3 }, 1000, _2nd, delta_f, 0, 1 }, 0.1 / 2, 5),
                          BiMaxPlasmaDesc({ { -O0, c / M_SQRT2, 2 }, 100, _1st, full_f }, 0.001 / 2));

    /// cold fluid plasma descriptors
    ///
    static constexpr auto cold_descs = std::make_tuple(ColdPlasmaDesc({ 0.000544662, 0.093352 }));

    /// external source descriptors
    ///
    static constexpr auto source_descs = std::make_tuple();

    //
    // MARK: Data Recording
    //

    /// a top-level directory to which outputs will be saved
    ///
    /// This setting is configurable through the command line option, `--wd <path_to_data_dump>`
    ///
    static constexpr char working_directory[] = "./data";

    /// field and particle energy density recording frequency; in units of inner_Nt
    /// `0' means `not interested'
    ///
    static constexpr unsigned energy_recording_frequency = 1;

    /// electric and magnetic field recording frequency
    ///
    static constexpr unsigned field_recording_frequency = 10;

    /// species moment recording frequency
    ///
    static constexpr unsigned moment_recording_frequency = 10;

    /// simulation particle recording frequency
    ///
    static constexpr unsigned particle_recording_frequency = 1000;

    /// maximum number of particles to dump
    ///
    static constexpr std::array<unsigned, std::tuple_size_v<decltype(part_descs)>> Ndumps
        = { 1000000, 1000000 };

    /// velocity histogram recording frequency
    ///
    static constexpr unsigned vhistogram_recording_frequency = 1000;

    /// per-species gyro-averaged velocity space specification used for sampling velocity histogram
    ///
    /// the parallel (v1) and perpendicular (v2) velocity specs are described by
    /// the range of the velocity space extent and the number of velocity bins
    ///
    /// note that the Range type is initialized with the OFFSET (or location) and the LENGTH
    ///
    /// recording histograms corresponding to specifications with the bin count being 0 will be
    /// skipped over
    ///
    static constexpr std::array<std::pair<Range, unsigned>, std::tuple_size_v<decltype(part_descs)>>
        v1hist_specs = { std::make_pair(1.3 * Range{ -1, 2 }, 150), std::make_pair(0.13 * Range{ -1, 2 }, 150) };
    static constexpr std::array<std::pair<Range, unsigned>, std::tuple_size_v<decltype(part_descs)>>
        v2hist_specs = { std::make_pair(1.3 * 2.3 * Range{ +0, 1 }, 80), std::make_pair(0.13 * Range{ +0, 1 }, 80) };
};

/// debugging options
///
namespace Debug {
constexpr bool zero_out_electromagnetic_field = false;
constexpr Real initial_efield_noise_amplitude = 0e0;
} // namespace Debug
