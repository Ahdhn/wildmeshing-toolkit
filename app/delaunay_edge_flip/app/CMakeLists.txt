# Include more application dependencies here
include(cli11)

file(GLOB APP_SOURCES CONFIGURE_DEPENDS "*.h" "*.cpp")
add_executable(def_app ${APP_SOURCES})
# Direct dependencies
target_link_libraries(def_app PUBLIC
    wmtk::toolkit
    wmtk::data
    igl::predicates
#        FastEnvelope
)
target_link_libraries(def_app PUBLIC wmtk_def_lib)
wmtk_copy_dll(def_app)

set_target_properties(def_app PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/app")
