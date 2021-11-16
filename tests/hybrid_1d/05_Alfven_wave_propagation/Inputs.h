/*
 * Copyright (c) 2021, Kyungguk Min
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
    static constexpr unsigned number_of_subdomains = 3;

    /// number of subdomain clones on which evenly divided particles are assigned and updated (positive integer)
    ///
    static constexpr unsigned number_of_distributed_particle_subdomain_clones = 2;

    /// number of worker threads to spawn for parallelization
    ///
    /// value `0' means serial update; value `n' means parallelization using n + 1 threads
    /// n + 1 must be divisible by number_of_subdomains * number_of_distributed_particle_subdomain_clones
    ///
    static constexpr unsigned number_of_worker_threads = 2 * number_of_subdomains * number_of_distributed_particle_subdomain_clones - 1;

    /// electric field extrapolation method
    ///
    static constexpr Algorithm algorithm = CAMCL;

    /// number of subscyles for magnetic field update; applied only for CAM-CL algorithm
    ///
    static constexpr unsigned n_subcycles = 10;

    /// particle boundary condition
    ///
    static constexpr BC particle_boundary_condition = BC::reflecting;

    /// wave masking function
    ///
    /// the first argument is masking inset, i.e., the number of grid points through which waves are gradually damped
    /// the second argument is the masking coefficients, zero being no masking at all and one being 0 to 100% masking within the masking inset
    ///
    /// see `docs/boundary_condition.nb` for how the phase retardation and amplitude damping work
    ///
    static constexpr MaskingFunction phase_retardation{ 100, 1.0 };
    static constexpr MaskingFunction amplitude_damping{ 100, 0.2 };

    //
    // MARK: Global parameters
    //

    /// light speed
    ///
    static constexpr Real c = 400;

    /// magnitude of equatorial background magnetic field
    ///
    static constexpr Real O0 = 1;

    /// inhomogeneity parameter, ξ
    /// the field variation at the central field line is given by B/B_eq = 1 + (ξ*x)^2,
    /// where x is the coordinate along the axis of mirror field symmetry
    ///
    static constexpr Real xi = 0.002;

    /// simulation grid size at the equator
    ///
    static constexpr Real Dx = 0.5;

    /// number of grid points
    ///
    static constexpr unsigned Nx = 1200;

    /// time step size
    ///
    static constexpr Real dt = 0.02;

    /// number of time steps for inner loop
    /// total time step Nt = inner_Nt * outer_Nt
    /// simulation time t = dt*Nt
    ///
    static constexpr unsigned inner_Nt = 50;

    /// number of time steps for outer loop
    /// total time step Nt = inner_Nt * outer_Nt
    /// simulation time t = dt*Nt
    ///
    /// This is configurable through the command line option, `--outer_Nt=[0-9].*`
    ///
    static constexpr unsigned outer_Nt = 1500;

    //
    // MARK: Plasma Species Descriptions
    //

    /// charge-neutralizing electron fluid description
    ///
    static constexpr auto efluid_desc = eFluidDesc({ -1836, 17139.4 });

    /// kinetic plasma descriptors
    ///
    static constexpr auto part_descs = std::make_tuple(
        // BiMaxPlasmaDesc({ { 1, 400, 1 }, 500, _2nd, delta_f, 0, 2 }, 0.001)
    );

    /// cold fluid plasma descriptors
    ///
    static constexpr auto cold_descs = std::make_tuple(ColdPlasmaDesc({ 1, 400 }));

    /// external source descriptors
    ///
    static constexpr auto source_descs = std::make_tuple(
        ExternalSourceDesc<1>({ 0.5, { 50, 2 * M_PI * 20 }, 30, 1 }, { ComplexVector{ 0, { 0, 0.04 * c }, { 0.04 * c, 0 } } }, { CurviCoord{ 100 } }));

    //
    // MARK: Data Recording
    //

    /// a top-level directory to which outputs will be saved
    ///
    static constexpr char working_directory[] = "./data";

    /// field and particle energy density recording frequency; in units of inner_Nt
    /// `0' means `not interested'
    ///
    static constexpr unsigned energy_recording_frequency = 1;

    /// electric and magnetic field recording frequency
    ///
    static constexpr unsigned field_recording_frequency = 1;

    /// ion species moment recording frequency
    ///
    static constexpr unsigned moment_recording_frequency = 0;

    /// simulation particle recording frequency
    ///
    static constexpr unsigned particle_recording_frequency = 0;

    /// maximum number of particles to dump
    ///
    static constexpr std::array<unsigned, std::tuple_size_v<decltype(part_descs)>> Ndumps
        = {};

    /// velocity histogram recording frequency
    ///
    static constexpr unsigned vhistogram_recording_frequency = 0;

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
        v1hist_specs = {};
    static constexpr std::array<std::pair<Range, unsigned>, std::tuple_size_v<decltype(part_descs)>>
        v2hist_specs = {};
};

/// debugging options
///
namespace Debug {
constexpr bool zero_out_electromagnetic_field = false;
constexpr Real initial_bfield_noise_amplitude = 0e0;
} // namespace Debug
