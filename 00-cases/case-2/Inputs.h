//
// code:    git@gitlab.com:space_plasma/mirror-pic_1d.git
// version: commit 0cc1812bcd55317a5ffeaa87c442509df2ada0cc
//

struct Input {
    //
    // MARK:- Environment
    //

    /// Number of ghost cells
    ///
    /// It must be greater than 0.
    ///
    static constexpr unsigned number_of_ghost_cells = 2;

    /// number of subdomains for domain decomposition (positive integer)
    ///
    /// Nx must be divisible by this number
    ///
    static constexpr unsigned number_of_subdomains = 7;

    /// number of subdomain clones on which evenly divided particles are assigned and updated (positive integer)
    ///
    static constexpr unsigned number_of_distributed_particle_subdomain_clones = 8 * 8;

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
    static constexpr BC particle_boundary_condition = BC::reflecting;

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
    static constexpr MaskingFunction phase_retardation{ 300, 1.0 };
    static constexpr MaskingFunction amplitude_damping{ 300, 0.2 };

    //
    // MARK: Global parameters
    //

    /// light speed
    ///
    static constexpr Real c = 1;

    /// magnitude of equatorial background magnetic field
    ///
    static constexpr Real O0 = 1;

    /// inhomogeneity parameter, ξ
    /// the field variation at the central field line is given by B/B_eq = 1 + (ξ*x)^2,
    /// where x is the coordinate along the axis of mirror field symmetry
    ///
    static constexpr Real xi = 0.00221359;

    /// simulation grid size at the equator
    ///
    static constexpr Real Dx = 0.07;

    /// number of grid points
    ///
    static constexpr unsigned Nx = 2800 * 2;

    /// time step size
    ///
    static constexpr Real dt = 0.01;

    /// number of time steps for inner loop
    /// total time step Nt = inner_Nt * outer_Nt
    /// simulation time t = dt*Nt
    ///
    static constexpr unsigned inner_Nt = 100;

    /// number of time steps for outer loop
    /// total time step Nt = inner_Nt * outer_Nt
    /// simulation time t = dt*Nt
    ///
    /// this option is configurable through the commandline option, e.g., "--outer_Nt=100"
    ///
    static constexpr unsigned outer_Nt = 10000;

    //
    // MARK: Plasma Species Descriptions
    //

    /// kinetic plasma descriptors
    ///
    static constexpr auto part_descs = std::make_tuple(
            LossconePlasmaDesc({ 0.3 }, { { -1, 0.316228, 3 }, 25000, _2nd, full_f }, 0.010125, 1.66123),
            PartialShellPlasmaDesc({ { -1, 0.316228, 3 }, 25000, _2nd, full_f }, 0.0002, 0, 0.16)
            );

    /// cold fluid plasma descriptors
    ///
    static constexpr auto cold_descs = std::make_tuple(
            ColdPlasmaDesc({ -1, 5 })
            );

    /// external source descriptors
    ///
    static constexpr Range ext_extent{ 100, 100000 };
    static constexpr auto source_descs = std::make_tuple(
            ExternalSourceDesc<1>({ 0.2, ext_extent, ext_extent.loc, 1 }, {
                ComplexVector{ 0, { 0, 0.008 * c }, { -0.008 * c, 0 } } }, { CurviCoord{ 0 } })
            );

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
    static constexpr std::pair<unsigned, Range> field_recording_frequency = { 1, { 100, 20000 } };

    /// a pair of
    ///
    /// - species moment recording frequency; in units of inner_Nt
    /// - recording start time and recording duration; in units of simulation time
    ///
    /// passing zero to the recording frequency means no recording
    ///
    /// the recording frequency option is configurable through the commandline option, e.g., "--moment_recording_frequency=10"
    ///
    static constexpr std::pair<unsigned, Range> moment_recording_frequency = { 1, { 100, 20000 } };

    /// a pair of
    ///
    /// - simulation particle recording frequency; in units of inner_Nt
    /// - recording start time and recording duration; in units of simulation time
    ///
    /// passing zero to the recording frequency means no recording
    ///
    /// the recording frequency option is configurable through the commandline option, e.g., "--particle_recording_frequency=10"
    ///
    static constexpr std::pair<unsigned, Range> particle_recording_frequency = { 50, { 3000, 6000 } };

    /// maximum number of particles to dump
    ///
    static constexpr std::array<unsigned, std::tuple_size_v<decltype(part_descs)>> Ndumps
        = { 70'000'000, 70'000'000 };

    /// a pair of
    ///
    /// - velocity histogram recording frequency; in units of inner_Nt
    /// - recording start time and recording duration; in units of simulation time
    ///
    /// passing zero to the recording frequency means no recording
    ///
    /// the recording frequency option is configurable through the commandline option, e.g., "--vhistogram_recording_frequency=10"
    ///
    static constexpr std::pair<unsigned, Range> vhistogram_recording_frequency = { 0, {} };

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
    static constexpr Real vA0 = 1.0 / 5.0;
    static constexpr std::array<std::pair<Range, unsigned>, std::tuple_size_v<decltype(part_descs)>>
        gv1hist_specs = {
            std::make_pair(Range{-1, 2} * 4.5 * vA0, 201),
            std::make_pair(Range{-1, 2} * 1.6 * vA0, 201),
        };
    static constexpr std::array<std::pair<Range, unsigned>, std::tuple_size_v<decltype(part_descs)>>
        gv2hist_specs = {
            std::make_pair(Range{+0, 1} * 5.5 * vA0, 101),
            std::make_pair(Range{+0, 1} * 1.6 * vA0, 101),
        };
};

/// debugging options
///
namespace Debug {
constexpr bool zero_out_electromagnetic_field = false;
constexpr Real initial_efield_noise_amplitude = 0e0;
constexpr bool should_use_unified_snapshot    = false;
} // namespace Debug
