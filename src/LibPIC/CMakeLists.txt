
add_library(LibPIC EXCLUDE_FROM_ALL)
#add_library(LibPIC)
target_sources(LibPIC PUBLIC
    ./PIC/Config.h
    ./PIC/Predefined.h
    ./PIC/println.h
    ./PIC/lippincott.h
    ./PIC/lippincott.cc
    ./PIC/Badge.h
    ./PIC/Options.h
    ./PIC/Options.cc
    ./PIC/Range.h
    ./PIC/Scalar.h
    ./PIC/Vector.h
    ./PIC/Tensor.h
    ./PIC/FourTensor.h
    ./PIC/FourVector.h
    ./PIC/Matrix.h
    ./PIC/CartCoord.h
    ./PIC/CurviCoord.h
    ./PIC/BorisPush.h
    ./PIC/MaskingFunction.h
    ./PIC/Shape.h
    ./PIC/Grid.h
    ./PIC/splitmix64.h
    ./PIC/xoroshiro128.h
    ./PIC/BitReversedPattern.h
    ./PIC/NRRandomEngine.h
    ./PIC/RandomReal.h
    ./PIC/RandomReal.cc
    ./PIC/MirrorBasis.h
    ./PIC/MirrorField.h
    ./PIC/MirrorTransform.h
    ./PIC/MirrorCotrans.h
    ./PIC/MirrorCotrans.cc
    ./PIC/MirrorGeometry.h
    ./PIC/MirrorGeometry.cc
    ./PIC/Geometry.h
    ./PIC/Particle.h
    ./PIC/RelativisticParticle.h
    ./PIC/TypeMaps.h
    ./PIC/PlasmaDesc.h
    ./PIC/VDF.h
    ./PIC/TestParticleVDF.h
    ./PIC/TestParticleVDF.cc
    ./PIC/MaxwellianVDF.h
    ./PIC/MaxwellianVDF.cc
    ./PIC/LossconeVDF.h
    ./PIC/LossconeVDF.cc
    ./PIC/PartialShellVDF.h
    ./PIC/PartialShellVDF.cc
    ./PIC/VDFVariant.h
    ./PIC/RelativisticVDF.h
    ./PIC/RelativisticTestParticleVDF.h
    ./PIC/RelativisticTestParticleVDF.cc
    ./PIC/RelativisticMaxwellianVDF.h
    ./PIC/RelativisticMaxwellianVDF.cc
    ./PIC/RelativisticLossconeVDF.h
    ./PIC/RelativisticLossconeVDF.cc
    ./PIC/RelativisticPartialShellVDF.h
    ./PIC/RelativisticPartialShellVDF.cc
    ./PIC/RelativisticVDFVariant.h
    )

set_target_properties(LibPIC PROPERTIES
    LINKER_LANGUAGE CXX
    OUTPUT_NAME PIC
    )
target_include_directories(LibPIC PUBLIC .)

set_project_warnings(LibPIC)
enable_sanitizers(LibPIC)

add_subdirectory(test)
