
add_library(LibPIC EXCLUDE_FROM_ALL)
#add_library(LibPIC)
target_sources(LibPIC PUBLIC
    PIC/RandomReal.cc
    PIC/lippincott.cc
    PIC/Options.cc
    PIC/MirrorGeometry.cc
    PIC/MirrorCotrans.cc
    PIC/TestParticleVDF.cc
    PIC/MaxwellianVDF.cc
    PIC/LossconeVDF.cc
    PIC/PartialShellVDF.cc
    PIC/RelativisticTestParticleVDF.cc
    PIC/RelativisticMaxwellianVDF.cc
    PIC/RelativisticLossconeVDF.cc
    PIC/RelativisticPartialShellVDF.cc
    )

set_target_properties(LibPIC PROPERTIES
    LINKER_LANGUAGE CXX
    OUTPUT_NAME PIC
    )
target_include_directories(LibPIC PUBLIC .)

set_project_warnings(LibPIC)
enable_sanitizers(LibPIC)

add_subdirectory(test)
