
Include(FetchContent)

FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v2.13.1)

FetchContent_MakeAvailable(Catch2)

add_executable(test-common)
target_sources(test-common PRIVATE
    main.cc
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
    )

target_link_libraries(test-common PRIVATE
    Catch2::Catch2 common
    )

set_project_warnings(test-common)
enable_sanitizers(test-common)

add_test("Test-libCommon" test-common)
