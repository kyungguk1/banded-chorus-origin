
#add_library(common EXCLUDE_FROM_ALL)
add_library(common)
target_sources(common PUBLIC
    common-config.h
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
#    VDF/VDF.h
#    VDF/MaxwellianVDF.h
#    VDF/LossconeVDF.h
#    TypeMaps.h
    )
target_sources(common PRIVATE
    Utility/lippincott.cc
    Utility/Options.cc
    Geometry.cc
#    VDF/VDF.cc
#    VDF/MaxwellianVDF.cc
#    VDF/LossconeVDF.cc
    )

set_target_properties(common PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(common PUBLIC .)

set_project_warnings(common)
enable_sanitizers(common)
