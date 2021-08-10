
#add_library(common EXCLUDE_FROM_ALL)
add_library(common)
target_sources(common PUBLIC
    common-config.h
    Utility/lippincott.h
    Utility/println.h
    Utility/Options.h
    Utility/Range.h
    Utility/Scalar.h
    Utility/Vector.h
    Utility/Tensor.h
    Utility/Shape.h
#    Utility/GridArray.h
    )
target_sources(common PRIVATE
    Utility/lippincott.cc
    Utility/Options.cc
#    Utility/BorisPush.h
#    Utility/Particle.h
#    Utility/TypeMaps.h
#    VDF/BitReversedPattern.h
#    VDF/BitReversedPattern.cc
#    VDF/MaxwellianVDF.h
#    VDF/MaxwellianVDF.cc
#    VDF/VDF.h
#    VDF/VDF.cc
#    VDF/LossconeVDF.h
#    VDF/LossconeVDF.cc
    )

set_target_properties(common PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(common PUBLIC .)

set_project_warnings(common)
enable_sanitizers(common)
