
add_library(PIC EXCLUDE_FROM_ALL)
#add_library(PIC)
target_sources(PIC PUBLIC
    PIC-config.h
    Predefined.h
    Utility/lippincott.h
    Utility/println.h
    Utility/Options.h
    Utility/Range.h
    Utility/Scalar.h
    Utility/Vector.h
    Utility/Tensor.h
    Utility/Shape.h
    Utility/Grid.h
    Utility/BorisPush.h
    Geometry.h
    PlasmaDesc.h
    VDF/Particle.h
    VDF/BitReversedPattern.h
    VDF/VDF.h
    VDF/MaxwellianVDF.h
    VDF/LossconeVDF.h
    VDF/VDFVariant.h
    Utility/TypeMaps.h
    )
target_sources(PIC PRIVATE
    Utility/lippincott.cc
    Utility/Options.cc
    Geometry.cc
    VDF/RandomReal.h
    VDF/MaxwellianVDF.cc
    VDF/LossconeVDF.cc
    )

set_target_properties(PIC PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(PIC PUBLIC .)

set_project_warnings(PIC)
enable_sanitizers(PIC)
