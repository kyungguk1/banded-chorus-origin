
add_library(common
    ./config.h
#    ./Utility/lippincott.h
#    ./Utility/lippincott.cc
#    ./Utility/Range.h
#    ./Utility/BorisPush.h
#    ./Utility/Shape.h
#    ./Utility/Scalar.h
#    ./Utility/Options.h
#    ./Utility/println.h
#    ./Utility/Particle.h
#    ./Utility/GridQ.h
#    ./Utility/Vector.h
#    ./Utility/Options.cc
#    ./Utility/Tensor.h
#    ./Utility/TypeMaps.h
#    ./VDF/BitReversedPattern.h
#    ./VDF/BitReversedPattern.cc
#    ./VDF/MaxwellianVDF.h
#    ./VDF/MaxwellianVDF.cc
#    ./VDF/VDF.h
#    ./VDF/VDF.cc
#    ./VDF/LossconeVDF.h
#    ./VDF/LossconeVDF.cc
    )

set_target_properties(common PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(common PUBLIC .)

set_project_warnings(common)
enable_sanitizers(common)
