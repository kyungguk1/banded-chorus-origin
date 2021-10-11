
add_library(PIC EXCLUDE_FROM_ALL)
#add_library(PIC)
target_sources(PIC PUBLIC
    PIC/Config.h
    PIC/Predefined.h
    PIC/Badge.h
    PIC/lippincott.h
    PIC/println.h
    PIC/Options.h
    PIC/Range.h
    PIC/Scalar.h
    PIC/Vector.h
    PIC/FourVector.h
    PIC/Tensor.h
    PIC/FourTensor.h
    PIC/Matrix.h
    PIC/CartCoord.h
    PIC/CurviCoord.h
    PIC/Shape.h
    PIC/Grid.h
    PIC/BorisPush.h
    PIC/PlasmaDesc.h
    PIC/Particle.h
    PIC/RelativisticParticle.h
    PIC/TypeMaps.h
    PIC/Geometry.h
    PIC/MirrorGeometry.h
    PIC/MirrorCotrans.h
    PIC/MirrorField.h
    PIC/MirrorBasis.h
    PIC/MirrorTransform.h
    PIC/VDF.h
    PIC/TestParticleVDF.h
    PIC/MaxwellianVDF.h
    PIC/LossconeVDF.h
    PIC/VDFVariant.h
    PIC/RelativisticVDF.h
    PIC/RelativisticTestParticleVDF.h
    PIC/RelativisticMaxwellianVDF.h
    PIC/RelativisticLossconeVDF.h
    PIC/RelativisticVDFVariant.h
    )
target_sources(PIC PRIVATE
    PIC/BitReversedPattern.h
    PIC/RandomReal.h
    PIC/splitmix64.h
    PIC/xoroshiro128.h
    PIC/lippincott.cc
    PIC/Options.cc
    PIC/MirrorGeometry.cc
    PIC/MirrorCotrans.cc
    PIC/TestParticleVDF.cc
    PIC/MaxwellianVDF.cc
    PIC/LossconeVDF.cc
    PIC/RelativisticTestParticleVDF.cc
    PIC/RelativisticMaxwellianVDF.cc
    PIC/RelativisticLossconeVDF.cc
    )

set_target_properties(PIC PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(PIC PUBLIC .)

set_project_warnings(PIC)
enable_sanitizers(PIC)

add_subdirectory(test)
