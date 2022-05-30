
add_executable(test-pic_1d EXCLUDE_FROM_ALL)
target_sources(test-pic_1d PRIVATE
    test-main.cc
    ../ParamSet.cc
    ../Core/Species.cc
    ../Core/ExternalSource.cc
    test-ParamSet.cc
    test-ExternalSource.cc
    )
target_include_directories(test-pic_1d PRIVATE ..)

target_link_libraries(test-pic_1d PRIVATE
    Catch2::Catch2 LibPIC
    )

set_project_warnings(test-pic_1d)
enable_sanitizers(test-pic_1d)

add_test("Test-pic_1d" test-pic_1d)
