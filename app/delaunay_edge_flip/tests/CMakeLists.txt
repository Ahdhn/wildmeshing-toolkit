include(catch2)
FetchContent_GetProperties(catch2)
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/contrib)
include(Catch)

file(GLOB TEST_SOURCES CONFIGURE_DEPENDS "*.h" "*.cpp")
add_executable(def_tests ${TEST_SOURCES})
target_link_libraries(def_tests PUBLIC
        wmtk::def_lib
        Catch2::Catch2
        )
wmtk_copy_dll(def_tests)

set_target_properties(def_tests PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests")

catch_discover_tests(def_tests)
