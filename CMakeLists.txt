cmake_minimum_required(VERSION 3.18)

project(mirror-pic_1d CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

include(cmake/StandardProjectSettings.cmake)
include(cmake/PreventInSourceBuilds.cmake)
include(cmake/CompilerWarnings.cmake)
include(cmake/StaticAnalyzers.cmake)
include(cmake/Sanitizers.cmake)
include(cmake/Doxygen.cmake)
set(WARNINGS_AS_ERRORS $ENV{WARNINGS_AS_ERRORS})

enable_doxygen()
enable_testing()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")
if (ENABLE_KMIN_BUILD_CONFIG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libc++")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -F/Library/com.kyungguk.Frameworks.localized")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -F/Library/com.kyungguk.Frameworks.localized")
    set(KMIN_MPILIB "-framework mpich2")
    set(KMIN_HDF5LIB "-framework hdf5")
else ()
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath,$ENV{LD_LIBRARY_PATH}")
    set(KMIN_HDF5LIB "-lhdf5")
endif ()

# fetch external libraries
Include(FetchContent)

FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v2.13.7
)
FetchContent_Declare(
    ParallelKit
    GIT_REPOSITORY https://gitlab.com/space_plasma/lib/ParallelKit.git
    GIT_TAG 573d92dc72f4ed5a937078930acae4ad939e6725
    SOURCE_SUBDIR ParallelKit # New in version 3.18
)
FetchContent_Declare(
    HDF5Kit
    GIT_REPOSITORY https://gitlab.com/space_plasma/lib/HDF5Kit.git
    GIT_TAG dda1f938d6fed074bbe34dfbcea0a9357ca6eaab
    SOURCE_SUBDIR HDF5Kit # New in version 3.18
)

FetchContent_MakeAvailable(Catch2 ParallelKit HDF5Kit)

target_link_libraries(ParallelKit
    PUBLIC "-lpthread"
    PUBLIC "${KMIN_MPILIB}"
    )
target_link_libraries(HDF5Kit
    PUBLIC "${KMIN_HDF5LIB}"
    )

# libPIC
add_subdirectory(src/libPIC)
target_link_libraries(PIC
    PUBLIC ParallelKit
    PUBLIC HDF5Kit
    )

# non-relativistic pic_1d
add_subdirectory(src/pic_1d)
target_link_libraries(pic_1d PRIVATE PIC)
if (DEFINED PIC_INPUT_DIR)
    target_include_directories(pic_1d BEFORE PRIVATE ${PIC_INPUT_DIR})
else ()
    target_include_directories(pic_1d PRIVATE src/pic_1d)
endif ()

# relativistic pic_1d
add_subdirectory(src/rel_pic_1d)
target_link_libraries(rel_pic_1d PRIVATE PIC)
target_include_directories(rel_pic_1d PRIVATE src/rel_pic_1d)
if (DEFINED PIC_INPUT_DIR)
    target_include_directories(rel_pic_1d BEFORE PRIVATE ${PIC_INPUT_DIR})
else ()
    target_include_directories(rel_pic_1d PRIVATE src/rel_pic_1d)
endif ()

# hybrid pic_1d
add_subdirectory(src/hybrid_1d)
target_link_libraries(hybrid_1d PRIVATE PIC)
target_include_directories(hybrid_1d PRIVATE src/hybrid_1d)
if (DEFINED HYBRID_INPUT_DIR)
    target_include_directories(hybrid_1d BEFORE PRIVATE ${HYBRID_INPUT_DIR})
else ()
    target_include_directories(hybrid_1d PRIVATE src/hybrid_1d)
endif ()
