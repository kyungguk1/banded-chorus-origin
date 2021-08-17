
add_library(PIC EXCLUDE_FROM_ALL)
#add_library(PIC)
target_sources(PIC PUBLIC
    PIC/Config.h
    PIC/Predefined.h
    PIC/lippincott.h
    PIC/println.h
    PIC/Options.h
    PIC/Range.h
    PIC/Scalar.h
    PIC/Vector.h
    PIC/Tensor.h
    PIC/Shape.h
    PIC/Grid.h
    PIC/BorisPush.h
    PIC/Geometry.h
    PIC/PlasmaDesc.h
    PIC/Particle.h
    PIC/BitReversedPattern.h
    PIC/VDF.h
    PIC/MaxwellianVDF.h
    PIC/LossconeVDF.h
    PIC/VDFVariant.h
    PIC/TypeMaps.h
    )
target_sources(PIC PRIVATE
    PIC/lippincott.cc
    PIC/Options.cc
    PIC/Geometry.cc
    PIC/RandomReal.h
    PIC/splitmix64.h
    PIC/xoroshiro128.h
    PIC/MaxwellianVDF.cc
    PIC/LossconeVDF.cc
    )

set_target_properties(PIC PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(PIC PUBLIC .)

set_project_warnings(PIC)
enable_sanitizers(PIC)

add_subdirectory(test)
