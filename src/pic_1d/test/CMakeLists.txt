
Include(FetchContent)

FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v2.13.1)

FetchContent_MakeAvailable(Catch2)

add_executable(test-pic_1d EXCLUDE_FROM_ALL)
target_sources(test-pic_1d PRIVATE
    test-main.cc
    ../ParamSet.h
    ../ParamSet.cc
    ../Core/Species.h
    ../Core/Species.cc
    ../Core/ExternalSource.h
    ../Core/ExternalSource.cc
    test-ParamSet.cc
    test-ExternalSource.cc
    )
target_include_directories(test-pic_1d PRIVATE ..)

target_link_libraries(test-pic_1d PRIVATE
    Catch2::Catch2 PIC
    )

set_project_warnings(test-pic_1d)
enable_sanitizers(test-pic_1d)

add_test("Test-pic_1d" test-pic_1d)
