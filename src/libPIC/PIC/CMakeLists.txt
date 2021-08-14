
add_library(PIC EXCLUDE_FROM_ALL)
#add_library(PIC)
target_sources(PIC PUBLIC
    Config.h
    Predefined.h
    lippincott.h
    println.h
    Options.h
    Range.h
    Scalar.h
    Vector.h
    Tensor.h
    Shape.h
    Grid.h
    BorisPush.h
    Geometry.h
    PlasmaDesc.h
    Particle.h
    BitReversedPattern.h
    VDF.h
    MaxwellianVDF.h
    LossconeVDF.h
    VDFVariant.h
    TypeMaps.h
    )
target_sources(PIC PRIVATE
    lippincott.cc
    Options.cc
    Geometry.cc
    RandomReal.h
    MaxwellianVDF.cc
    LossconeVDF.cc
    )

set_target_properties(PIC PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(PIC PUBLIC ..)

set_project_warnings(PIC)
enable_sanitizers(PIC)
