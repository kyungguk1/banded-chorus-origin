/*
 * Copyright (c) 2019-2022, Kyungguk Min
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

    /// Number of ghost cells
    ///
    /// It must be greater than 0.
    ///
    static constexpr unsigned number_of_ghost_cells = 3;

    /// number of subdomains for domain decomposition (positive integer)
    ///
    /// Nx must be divisible by this number
    ///
    static constexpr unsigned number_of_subdomains = 4;

    /// number of subdomain clones on which evenly divided particles are assigned and updated (positive integer)
    ///
    static constexpr unsigned number_of_distributed_particle_subdomain_clones = 3;

    /// number of worker threads to spawn for parallelization
    ///
    /// value `0' means serial update; value `n' means parallelization using n + 1 threads
    /// n + 1 must be divisible by number_of_subdomains * number_of_distributed_particle_subdomain_clones
    ///
    static constexpr unsigned number_of_worker_threads
        = 1 * number_of_subdomains * number_of_distributed_particle_subdomain_clones - 1;

    /// flag to suppress longitudinal and/or transverse components of the field fluctuations
    ///
    /// setting both to true will zero out all field fluctuations
    ///
    static constexpr bool should_neglect_transverse_component   = false;
    static constexpr bool should_neglect_longitudinal_component = true;

    /// particle boundary condition
    ///
    static constexpr BC particle_boundary_condition = BC::periodic;

    /// if set, randomize the gyro-phase of reflected particles
    ///
    static constexpr bool should_randomize_gyrophase_of_reflecting_particles = true;

    /// wave masking function
    ///
    /// the first argument is masking inset, i.e., the number of grid points through which waves are gradually damped
    /// the second argument is the masking coefficients, zero being no masking at all and one being 0 to 100% masking within the masking inset
    ///
    /// see `docs/boundary_condition.nb` for how the phase retardation and amplitude damping work
    ///
    static constexpr MaskingFunction phase_retardation{};
    static constexpr MaskingFunction amplitude_damping{};

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
    static constexpr Real xi = 0;

    /// simulation grid size at the equator
    ///
    static constexpr Real Dx = 0.3;

    /// number of grid points
    ///
    static constexpr unsigned Nx = 480;

    /// time step size
    ///
    static constexpr Real dt = 0.02;

    /// number of time steps for inner loop
    /// total time step Nt = inner_Nt * outer_Nt
    /// simulation time t = dt*Nt
    ///
    static constexpr unsigned inner_Nt = 1;

    /// number of time steps for outer loop
    /// total time step Nt = inner_Nt * outer_Nt
    /// simulation time t = dt*Nt
    ///
    /// this option is configurable through the commandline option, e.g., "--outer_Nt=100"
    ///
    static constexpr unsigned outer_Nt = 1;

    //
    // MARK: Plasma Species Descriptions
    //

    /// kinetic plasma descriptors
    ///
    static constexpr auto part_descs = std::make_tuple(
        BiMaxPlasmaDesc({ { -1, 4, 2 }, 10000, TSC, full_f }, 1, 2),
        BiMaxPlasmaDesc({ { 0.000544662, 0.093352, 2 }, 10000, TSC, full_f }, 0.01));

    /// cold fluid plasma descriptors
    ///
    static constexpr auto cold_descs = std::make_tuple();

    /// external source descriptors
    ///
    static constexpr auto source_descs = std::make_tuple(
        ExternalSourceDesc<1>{ { 1, { 10, 100 }, 3, 2 }, { ComplexVector{ 1, 1i, 0 } }, { CurviCoord{ Nx + 1 } } });

    //
    // MARK: Data Recording
    //

    /// a top-level directory to which outputs will be saved
    ///
    /// this option is configurable through the commandline option, e.g., "--wd ./data"
    ///
    static constexpr std::string_view working_directory = "./data";

    /// a pair of
    ///
    /// - field and particle energy density recording frequency; in units of inner_Nt
    /// - recording start time and recording duration; in units of simulation time
    ///
    /// passing zero to the recording frequency means no recording
    ///
    /// the recording frequency option is configurable through the commandline option, e.g., "--energy_recording_frequency=10"
    ///
    static constexpr std::pair<unsigned, Range> energy_recording_frequency = { 1, {} };

    /// a pair of
    ///
    /// - electric and magnetic field recording frequency; in units of inner_Nt
    /// - recording start time and recording duration; in units of simulation time
    ///
    /// passing zero to the recording frequency means no recording
    ///
    /// the recording frequency option is configurable through the commandline option, e.g., "--field_recording_frequency=10"
    ///
    static constexpr std::pair<unsigned, Range> field_recording_frequency = { 5, {} };

    /// a pair of
    ///
    /// - species moment recording frequency; in units of inner_Nt
    /// - recording start time and recording duration; in units of simulation time
    ///
    /// passing zero to the recording frequency means no recording
    ///
    /// the recording frequency option is configurable through the commandline option, e.g., "--moment_recording_frequency=10"
    ///
    static constexpr std::pair<unsigned, Range> moment_recording_frequency = { 5, {} };

    /// a pair of
    ///
    /// - simulation particle recording frequency; in units of inner_Nt
    /// - recording start time and recording duration; in units of simulation time
    ///
    /// passing zero to the recording frequency means no recording
    ///
    /// the recording frequency option is configurable through the commandline option, e.g., "--particle_recording_frequency=10"
    ///
    static constexpr std::pair<unsigned, Range> particle_recording_frequency = { 1, {} };

    /// maximum number of particles to dump
    ///
    static constexpr std::array<unsigned, std::tuple_size_v<decltype(part_descs)>> Ndumps
        = { 10000 };

    /// a pair of
    ///
    /// - velocity histogram recording frequency; in units of inner_Nt
    /// - recording start time and recording duration; in units of simulation time
    ///
    /// passing zero to the recording frequency means no recording
    ///
    /// the recording frequency option is configurable through the commandline option, e.g., "--vhistogram_recording_frequency=10"
    ///
    static constexpr std::pair<unsigned, Range> vhistogram_recording_frequency = { 1, {} };

    /// per-species gyro-averaged momentum space specification used for sampling momentum histogram
    ///
    /// the parallel (γ*v1) and perpendicular (γ*v2) momentum specs are described by
    /// the range of the momentum space extent and the number of momentum bins
    ///
    /// note that the Range type is initialized with an OFFSET (or location) and LENGTH
    ///
    /// recording histograms corresponding to specifications with the bin count being 0 will be
    /// skipped over
    ///
    static constexpr std::array<std::pair<Range, unsigned>, std::tuple_size_v<decltype(part_descs)>>
        gv1hist_specs = { std::make_pair(Range{ -3, 6 }, 200) };
    static constexpr std::array<std::pair<Range, unsigned>, std::tuple_size_v<decltype(part_descs)>>
        gv2hist_specs = { std::make_pair(Range{ 0, 4 }, 100) };
};

/// debugging options
///
namespace Debug {
constexpr bool zero_out_electromagnetic_field = false;
constexpr Real initial_efield_noise_amplitude = 0e0;
constexpr bool should_use_unified_snapshot    = false;
} // namespace Debug
