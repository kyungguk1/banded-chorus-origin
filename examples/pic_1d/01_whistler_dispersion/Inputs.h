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
    // MARK:- Housekeeping
    //

    /// number of subdomains for domain decomposition (positive integer)
    ///
    /// Nx must be divisible by this number
    ///
    static constexpr unsigned number_of_subdomains = 10;

    /// number of worker threads to spawn for parallelization
    ///
    /// value `0' means serial update; value `n' means parallelization using n + 1 threads
    /// part_desc.Nc*Nx must be divisible by n + 1, and
    /// n + 1 must be divisible by number_of_subdomains
    ///
    static constexpr unsigned number_of_worker_threads = number_of_subdomains - 1;

    /// flag to suppress transverse electromagnetic fields
    ///
    static constexpr bool is_electrostatic = false;

    //
    // MARK: Global parameters
    //

    /// light speed
    ///
    static constexpr Real c = 4;

    /// magnitude of uniform background magnetic field
    ///
    static constexpr Real O0 = 1;

    /// angle in degrees between the x-axis and the uniform magnetic field direction
    ///
    static constexpr Real theta = 0;

    /// simulation grid size
    ///
    static constexpr Real Dx = 0.3 / 5;

    /// number of grid points
    ///
    static constexpr unsigned Nx = 960 * 2;

    /// time step size
    ///
    static constexpr Real dt = 0.02 / 5;

    /// number of time steps for inner loop
    ///
    /// total time step Nt = inner_Nt * outer_Nt
    /// simulation time t = dt*Nt
    ///
    static constexpr unsigned inner_Nt = 15 * 2;

    /// number of time steps for outer loop
    /// total time step Nt = inner_Nt * outer_Nt
    /// simulation time t = dt*Nt
    ///
    /// This is configurable through the command line option, `--outer_Nt=[0-9].*`
    ///
    static constexpr unsigned outer_Nt = 960;

    //
    // MARK: Plasma Species Descriptions
    //

    /// kinetic plasma descriptors
    ///
    static constexpr auto part_descs
        = std::make_tuple(BiMaxPlasmaDesc({ { -1, 4, 2 }, 1000, TSC, full_f }, 0.01));

    /// cold fluid plasma descriptors
    ///
    static constexpr auto cold_descs = std::make_tuple(ColdPlasmaDesc({ 0.000544662, 0.093352 }));

    //
    // MARK: Data Recording
    //

    /// a top-level directory to which outputs will be saved
    ///
    /// This setting is configurable through the command line option, `--wd <path_to_data_dump>`
    ///
    static constexpr char working_directory[] = "./data";

    /// field and particle energy density recording frequency; in units of inner_Nt
    ///
    static constexpr unsigned energy_recording_frequency = 1;

    /// electric and magnetic field recording frequency
    ///
    static constexpr unsigned field_recording_frequency = 1;

    /// species moment recording frequency
    ///
    static constexpr unsigned moment_recording_frequency = 1;

    /// simulation particle recording frequency
    ///
    static constexpr unsigned particle_recording_frequency = 1000;

    /// maximum number of particles to dump
    ///
    static constexpr std::array<unsigned, std::tuple_size_v<decltype(part_descs)>> Ndumps
        = { ~(0U) };

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
        v1hist_specs = { std::make_pair(0.4 * Range{ -1, 2 }, 150) };
    static constexpr std::array<std::pair<Range, unsigned>, std::tuple_size_v<decltype(part_descs)>>
        v2hist_specs = { std::make_pair(0.4 * Range{ +0, 1 }, 80) };
};

/// debugging options
///
namespace Debug {
constexpr bool zero_out_electromagnetic_field = false;
constexpr Real initial_efield_noise_amplitude = 0e0;
} // namespace Debug
