
add_library(PIC EXCLUDE_FROM_ALL)
#add_library(PIC)
target_sources(PIC PUBLIC
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

set_target_properties(PIC PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(PIC PUBLIC .)

set_project_warnings(PIC)
enable_sanitizers(PIC)

add_subdirectory(test)
