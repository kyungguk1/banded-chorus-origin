# Ion Cyclotron Wave Dispersion Relation

This example shows the ion cyclotron mode dispersion relation in a plasma with dominant protons and minor helium ions.

## Build Instruction

Make sure you have a working C++ compiler and a MPI library.

1. Create a build directory:

```
mkdir build && cd build
```

2. Configure CMake build system:

```
cmake -DHYBRID_INPUT_DIR=$PWD/.. \
-DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_BUILD_TYPE=Release \
-DENABLE_IPO=On -G "Ninja" ../../../..
```

The `HYBRID_INPUT_DIR` variable specifies the directory containing the `Inputs.h` header, which is located in the parent
directory. The last argument `../../../..` is the project directory where the top-level `CMakeLists.txt` is located. If
you don't have `ninja` installed but have `make`, replace `"Ninja"` with `"UNIX Makefiles"`.

3. Build:

```
ninja hybrid_1d
```

If you have used `"UNIX Makefiles"` in Step 3, replace `ninja` with `make`. `ninja` automatically parallelize the source
code compilation. You can specify the number of parallelism with an option, e.g., `-j4`, which instructs four
simultaneous compilations. The `j` option also applies to the `make` command.

## Simulation Configuration

> Whenever you change the settings here, you need to rebuild the program in Step 3 before execution.

### Global Configuration

All the input parameters are contained in the `Inputs.h` header. In this example, the speed is normalized to the proton
Alfven speed (assuming all ions are protons); the time to the inverse of the proton cyclotron frequency; and the length
to the proton inertial length.

* The `number_of_subdomains` parameter specifies the number of MPI processes to spawn. The total number of grid points
  must be divisible by this number. Internally, the simulation domain is broken up into this many *subdomains*. One MPI
  process is responsible for one subdomain.
* The `algorithm` parameter specifies the algorithm to extrapolate the electric field. Supported options are `CAMCL` (
  Current Advance Method with Cyclic Leapfrog) and `PC` (Predictor-Corrector).
* The `n_subcycles` parameter specifies the number of sub-cycles to update the magnetic field. This only applies to
  the `CAMCL` algorithm.
* The `c` parameter specifies the light speed. (The light propagation does not come into play in a hybrid code; this is
  merely for normalization.)
* The `O0` parameter specifies the magnitude of the uniform background magnetic field. The magnetic and electric fields
  of the simulation outputs are normalized to this value.
* The `theta` parameter specifies the angle, in degrees, between the *x*-axis and the background magnetic field vector.
* The `Dx` parameter specifies the simulation grid size, i.e., spatial resolution.
* The `Nx` parameter specifies the number of grid points. The whole simulation domain length is then `Lx = Dx * Nx`.
* The `dt` parameter specifies the simulation time step size.
* The `inner_Nt` parameter specifies the inner for-loop count.
* The `outer_Nt` parameter specifies the outer for-loop count. This parameter is configurable through the command line
  option, `--outer_Nt=<non-negative-integer>`.

> The system is evolved in time in a nested for-loop.
> The inner loop iterates `inner_Nt` times before exiting to the outer loop, where
> the program does some housekeeping tasks, such as recording the global iteration count
> and saving intermediate outputs (if there are any). Once all is done, it enters the inner loop and
> the cycle continues until the global iteration count reaches `inner_Nt * outer_Nt`.

### Plasma Configuration

The `efluid_desc` parameter specifies the charge-neutralizing, massless electron fluid. In this example,

```
eFluidDesc({ -1836, 9180.01 }, 0.01, isothermal)
```

The first parameter inside the curly braces specifies the electron cyclotron frequency (which carries the sign of the
electron charge). The second parameter inside the curly braces specifies the electron plasma frequency. The second
parameter to `eFluidDesc` species the electron beta. The last parameter specifies the polytropic index:
either `isothermal` or `adiabatic`.

The `part_descs` parameter contains information about kinetic plasma populations. In this example, There is two
populations: dense protons and tenuous helium ions. Both populations are Maxwellian-distributed. All information needed
to construct a Maxwellian is passed by an object, `BiMaxPlasmaDesc`. In this example,

```
BiMaxPlasmaDesc({ { 1, 213.169, 2 }, 1000, TSC, full_f }, 0.0099, 1, 0),
BiMaxPlasmaDesc({ { .25, 10.7122, 2 }, 1000, TSC, full_f }, 0.0001, 1, 0)
```

The first parameter in the inner curly braces specifies the cyclotron frequency. Since a helium ion is four times
heavier than a proton, the cyclotron frequency of the second component is one quarter of the first one. The second
parameter specifies the plasma frequency. Note that this takes into account 99% concentration of protons and 1%
concentration of helium ions. For example, this parameter will be the same as `c` if all ion species are 100% protons.
Likewise, it will be 0.5*`c` if all ion species are helium ions. The last one, `2`, specifies the number of source (
i.e., charge density and current) smoothing; this is to reduce the discrete particle noise.

In the outer curly braces, the number `1000` specifies the number of simulation particles per grid cell. The `TSC`, a
synonym for `_2nd`, denotes the order of the particle shape function. The `full_f` parameter specifies to use the
*full-f* approach. The parameter next to the outer closing curly brace is plasma beta. The temperatures of both species
are assumed to be the same. The different number shown here is due to the relative concentrations of the two species.

The `cold_descs` parameter contains information about cold, linearized fluid populations. In this example, nothing is
specified.

### Intermediate Output Configuration

* The `working_directory` parameter specifies the relative, or absolute, path to a directory under which all simulation
  outputs should be saved. This parameter is configurable through the command line option, `--wd <path_to_save_dir>`.
* The `energy_recording_frequency` parameter specifies the frequency with which the average energy densities of various
  physical parameters should be saved. The value `1` means at every `inner_Nt` iterations. `0` means disabling output of
  this data product.
* The `field_recording_frequency` parameter specifies the frequency with which the electric and magnetic fields in the
  whole simulation domain should be saved.
* The `moment_recording_frequency` parameter specifies the frequency with which the particle (including cold fluid)
  moments in the whole simulation domain should be saved.
* The `particle_recording_frequency` parameter specifies the frequency with which a subset of the simulation particles
  should be saved.
* The `Ndumps` parameter is an array of the same length as the number of the kinetic plasma populations. Each element
  specifies the (rough) number of simulation particles out of the corresponding plasma population should be saved. The
  absence of numbers in this example means `0`.
* The `vhistogram_recording_frequency` parameter specifies the frequency with which the two-dimensional histogram of the
  simulation particles in the whole simulation domain in parallel (*v1*) and perpendicular (*v2*) velocity space should
  be saved.
* The `v1hist_specs` and `v2hist_specs` parameters specify the binning specification of a two-dimensional histogram.
  Each of them is an array of pairs of the same length. The length of the array must not exceed the number of kinetic
  plasma populations. Each pair is made of an extent of binning space (wrapped in a `Range` object) and the number of
  bins. In the example, the binning specification in *v1* space is

```
v1hist_specs = { std::make_pair(0.4 * Range{ -1, 2 }, 150) };
```

It is a one-element array, corresponding to the initially Maxwellian-distributed electron population. The `Range` object
expects the starting location as its first argument and the length of the range as its second argument. In addition,
simple arithmetic operations (`+`, `-`, and `*`) with a scalar on a *Range* object are allowed. So, the first element
in `v1hist_specs` means 150 bins evenly-spaced in the range *-0.4 < v1 < 0.4*.

Likewise, the binning specification in *v2* space

```
v2hist_specs = { std::make_pair(0.4 * Range{ +0, 1 }, 80) };
```

means 80 bins evenly-spaced in the range *0 < v1 < 0.4*.

## Simulation Run

In a terminal,

```
mpiexec -n 10 ./src/hybrid_1d/hybrid_1d --wd=<path_to_save_dir> --output_Nt=1000 --save=false --load=false
```

The number after `-n` option is the number of MPI processes to spawn, which must be the same as `number_of_subdomains`
in the global configuration. `./src/hybrid_1d/hybrid_1d` is the executable built. The `--save` and `--load` options are
optional, whose default values are `false`. When `--save=true` is set (a shorthand is `-save`), the program snapshots
the last state of the simulation at the end of the simulation run. When `--load` is set (a shorthand is `-load`), the
simulation resumes from the last saved state. Therefore, this assumes a snapshot already exists. A typical use of
the `-save` and `-load` combination is illustrated below:

```
$ mpiexec -n 10 ./src/hybrid_1d/hybrid_1d --wd ~/Downloads/test2 --outer_Nt 1 -save
Driver> initializing domain(s)
	initializing particles
master_loop> steps(x15) = 1/1; time = 0
%% time elapsed: 0.250887s
	saving snapshots
$ mpiexec -n 10 ./src/hybrid_1d/hybrid_1d --wd ~/Downloads/test2 --outer_Nt 100 -save -load
Driver> initializing domain(s)
	loading snapshots
master_loop> steps(x15) = 1/100; time = 0.3
master_loop> steps(x15) = 2/100; time = 0.6
master_loop> steps(x15) = 3/100; time = 0.9
...
master_loop> steps(x15) = 98/100; time = 29.4
master_loop> steps(x15) = 99/100; time = 29.7
master_loop> steps(x15) = 100/100; time = 30
%% time elapsed: 24.0811s
	saving snapshots
```

At first, the `--outer_Nt` option is set to 1 and `-save` is set. After program exits, a snapshot is saves and prints
out the time elapsed. You can use this time estimation to choose the number of iterations for a given time. In
subsequent runs, `-load` is also set to resume the simulation with a new number for the `--outer_Nt` option.
