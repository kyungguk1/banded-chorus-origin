add_executable(test-LibPIC EXCLUDE_FROM_ALL)
target_sources(test-LibPIC PRIVATE
    test-main.cc
    test-Badge.cc
    test-Options.cc
    test-Range.cc
    test-Scalar.cc
    test-Vector.cc
    test-Tensor.cc
    test-Matrix.cc
    test-Shape.cc
    test-Grid.cc
    test-BorisPush.cc
    test-Particle.cc
    test-RelativisticParticle.cc
    test-BitReversedPattern.cc
    test-Geometry.cc
    test-MirrorGeometry.cc
    test-PlasmaDesc.cc
    test-TypeMaps.cc
    test-xoroshiro128.cc
    test-NRRandomEngine.cc
    test-RandomReal.cc
    test-VDF.cc
    test-VDFVariant.cc
    test-RelativisticVDF.cc
    test-RelativisticVDFVariant.cc
    test-MaskingFunction.cc
    )

target_link_libraries(test-LibPIC PRIVATE
    Catch2::Catch2 LibPIC
    )

set_project_warnings(test-LibPIC)
enable_sanitizers(test-LibPIC)

add_test("Test-LibPIC" test-LibPIC)
