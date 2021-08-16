# Toy 1D Particle-in-cell Codes

## Preface

This project hosts two types of particle-in-cell plasma simulation codes popular in space physics community:

- Full particle-in-cell (PIC) code that treats both ions and electrons kinetically; and
- Hybrid code where ions are kinetic particles and electrons are a massless fluid.

In both cases, the spatial simulation domain is one-dimensional, but all three components of vector quantities are
retained. In addition, periodic boundary conditions are used.

The purpose of this project includes

- Test new idea
- Develop teaching materials
- Use as proving ground
- etc.

## Build Instruction

The project uses CMake as a build environment. All codes are written in C++. Building the targets requires MPI library
and hdf5 libraries.

Follow the steps below to build executables:

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
CXX=mpicxx cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_IPO=On \
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
