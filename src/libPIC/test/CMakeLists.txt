
Include(FetchContent)

FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v2.13.1)

FetchContent_MakeAvailable(Catch2)

add_executable(test-libPIC EXCLUDE_FROM_ALL)
target_sources(test-libPIC PRIVATE
    test-main.cc
    test-Options.cc
    test-Range.cc
    test-Scalar.cc
    test-Vector.cc
    test-Tensor.cc
    test-Shape.cc
    test-Grid.cc
    test-BorisPush.cc
    test-Particle.cc
    test-BitReversedPattern.cc
    test-Geometry.cc
    test-PlasmaDesc.cc
    test-VDF.cc
    test-VDFVariant.cc
    test-TypeMaps.cc
    test-RandomGenerator.cc
    ../PIC/splitmix64.h
    ../PIC/xoroshiro128.h
    )

target_link_libraries(test-libPIC PRIVATE
    Catch2::Catch2 PIC
    )

set_project_warnings(test-libPIC)
enable_sanitizers(test-libPIC)

add_test("Test-libPIC" test-libPIC)
