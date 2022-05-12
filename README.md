# Archive of the PIC code used for banded chorus simulations

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
    -G "Ninja" ${PROJECT_PATH}
```

- If `ninja` is not available, replace `"Ninja"` with `"Unix Makefiles"`.
- `PROJECT_PATH` refers to the project directory you just cloned.
- Set `PATH_TO_PIC_SIMULATION_INPUT_HEADER` to the path to a directory containing `Input.h`, if you are running the full
particle-in-cell code. Otherwise, exclude the whole `PIC_INPUT_DIR` option.
hyrid code. Otherwise, exclude the whole `HYBRID_INPUT_DIR` option.

4. Build the executables

```shell
ninja
```

The executable built is available at `src/rel_pic_1d/rel_pic_1d`.
