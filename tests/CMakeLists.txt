
Include(FetchContent)

FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v2.13.1)

FetchContent_MakeAvailable(Catch2)

add_executable(test-common EXCLUDE_FROM_ALL)
target_sources(test-common PRIVATE
    main.cc test-Options.cc
    )

target_link_libraries(test-common PRIVATE
    Catch2::Catch2 common
    )

set_project_warnings(test-common)
enable_sanitizers(test-common)

add_test(test-common test-common)
