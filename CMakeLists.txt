cmake_minimum_required(VERSION 3.17)

project(toy-pic_1d CXX)
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
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath,$ENV{LD_LIBRARY_PATH}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -F/Library/com.kyungguk.Frameworks.localized")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -F/Library/com.kyungguk.Frameworks.localized")
set(KMIN_MPILIB "-framework mpich2")
set(KMIN_HDF5LIB "-framework hdf5")

# external library
add_subdirectory(lib/ParallelKit/ParallelKit)
target_link_libraries(ParallelKit
    PUBLIC "-lpthread"
    PUBLIC "${KMIN_MPILIB}"
    )
add_subdirectory(lib/HDF5Kit/HDF5Kit)
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
if (DEFINED REL_PIC_INPUT_DIR)
    target_include_directories(rel_pic_1d BEFORE PRIVATE ${REL_PIC_INPUT_DIR})
else ()
    target_include_directories(rel_pic_1d PRIVATE src/rel_pic_1d)
endif ()

# hybrid pic_1d
add_subdirectory(src/hybrid_1d)
target_link_libraries(hybrid_1d
    PRIVATE ParallelKit
    PRIVATE HDF5Kit
    )
target_include_directories(hybrid_1d PRIVATE src/hybrid_1d)
if (DEFINED HYBRID_INPUT_DIR)
    target_include_directories(hybrid_1d BEFORE PRIVATE ${HYBRID_INPUT_DIR})
endif ()
