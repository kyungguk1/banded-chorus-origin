cmake_minimum_required(VERSION 3.17)

project(pic_1d CXX)
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
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libc++")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -F/Library/com.kyungguk.Frameworks.localized")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -F/Library/com.kyungguk.Frameworks.localized")
set(KMIN_MPILIB "-framework mpich2")
set(KMIN_HDF5LIB "-framework hdf5")

add_subdirectory(lib/ParallelKit/ParallelKit)
target_link_libraries(ParallelKit
    PUBLIC "-lpthread"
    )
add_subdirectory(lib/HDF5Kit/HDF5Kit)

add_subdirectory(src/pic_1d)
target_link_libraries(pic_1d
    PRIVATE ParallelKit
    PRIVATE "${KMIN_MPILIB}"
    PRIVATE HDF5Kit
    PRIVATE "${KMIN_HDF5LIB}"
    )
target_include_directories(pic_1d PRIVATE src/pic_1d)
if (DEFINED PIC_INPUT_DIR)
    target_include_directories(pic_1d BEFORE PRIVATE ${PIC_INPUT_DIR})
endif ()

add_subdirectory(src/hybrid_1d)
target_link_libraries(hybrid_1d
    PRIVATE ParallelKit
    PRIVATE "${KMIN_MPILIB}"
    PRIVATE HDF5Kit
    PRIVATE "${KMIN_HDF5LIB}"
    )
target_include_directories(hybrid_1d PRIVATE src/hybrid_1d)
if (DEFINED HYBRID_INPUT_DIR)
    target_include_directories(hybrid_1d BEFORE PRIVATE ${HYBRID_INPUT_DIR})
endif ()
