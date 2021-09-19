# 1D Particle-in-cell Code in Mirror Magnetic Field

This project hosts two types of particle-in-cell (PIC) codes written in C++:

- Full particle-in-cell code that treats both ions and electrons kinetically; and
- Hybrid code where ions are kinetic particles and electrons are a massless fluid.

Both are one-dimensional along a mirror background magnetic field.

## Build Instruction

The project uses CMake as a build environment. Building the targets requires MPI and hdf5 as an external dependence.
Avoid compiling with OpenMPI.

Follow the steps below to build executables:

0. Clone the project

1. Update submodules inside the cloned project root directory

```
git submodule update --init
```

2. Make a build directory in the directory where simulations are to be performed

```
mkdir build && cd build
```

3. Generate build configurations

```
cmake -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_BUILD_TYPE=Release -DENABLE_IPO=On \
-DPIC_INPUT_DIR=${PATH_TO_PIC_SIMULATION_INPUT_HEADER} \
-DHYBRID_INPUT_DIR=${PATH_TO_HYBRID_SIMULATION_INPUT_HEADER} \
-G "Ninja" ${PROJECT_PATH}
```

`PATH_TO_PIC_SIMULATION_INPUT_HEADER` is set to the path to a directory containing `Input.h` for PIC simulations; and
`PATH_TO_HYBRID_SIMULATION_INPUT_HEADER` is set to the path to a directory containing `Input.h` for hybrid simulations.

4. Build executables

```
ninja ${TARGET}
```

Here, `TARGET` is either `pic_1d`, `rel_pic_1d`, or `hybrid_1d`. In Step 3, the `HYBRID_INPUT_DIR` configuration option
can be omitted for the first two targets, while the `PIC_INPUT_DIR` configuration option can be omitted for the last
target.
