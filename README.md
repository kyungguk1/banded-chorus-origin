# 1D Particle-in-cell Code in Mirror Magnetic Field

This project hosts two types of particle-in-cell (PIC) codes written in C++:

- Full particle-in-cell code that treats both ions and electrons kinetically; and
- Hybrid code where ions are kinetic particles and electrons are a massless fluid.

Both are one-dimensional and take into account mirror-like magnetic field geometry.

## Directory Structure

- The `docs` directory contains documents of the numerical implementation and mathematical derivations.
- The `tests` directory contains simple example simulation setups.
- The `src` directory contains the source codes.

## Build Instruction

The project uses CMake as a build environment. Building the targets requires **MPI** and **hdf5** libraries available.

> Avoid compiling with OpenMPI. This implementation seems to misbehave when communicating struct-like data.

Once all the dependencies are met, follow the steps below to build executables:

1. Clone the project

Follow the instruction in the project page.

2. Make a build directory

```shell
mkdir build && cd build
```

3. Generate the build configurations

```shell
cmake -DCMAKE_CXX_COMPILER=mpicxx \
    -DCMAKE_BUILD_TYPE=Release -DENABLE_IPO=On \
    -DPIC_INPUT_DIR=${PATH_TO_PIC_SIMULATION_INPUT_HEADER} \
    -DHYBRID_INPUT_DIR=${PATH_TO_HYBRID_SIMULATION_INPUT_HEADER} \
    -G "Ninja" ${PROJECT_PATH}
```

- If `ninja` is not available, replace `"Ninja"` with `"Unix Makefiles"`.
- `PROJECT_PATH` refers to the project directory you just cloned.
- Set `PATH_TO_PIC_SIMULATION_INPUT_HEADER` to the path to a directory containing `Input.h`, if you are running the full
particle-in-cell code. Otherwise, exclude the whole `PIC_INPUT_DIR` option.
- Set `PATH_TO_HYBRID_SIMULATION_INPUT_HEADER` to the path to a directory containing `Input.h`, if you are running the
hyrid code. Otherwise, exclude the whole `HYBRID_INPUT_DIR` option.

4. Build the executables

```shell
ninja ${TARGET}
```

- `TARGET` is either `pic_1d`, `rel_pic_1d`, or `hybrid_1d`.
- If `"Unix Makefiles"` has been used in the configuration phase, replace `ninja` with `make`.

The executable built is available at `src/${TARGET}/${TARGET}`.
