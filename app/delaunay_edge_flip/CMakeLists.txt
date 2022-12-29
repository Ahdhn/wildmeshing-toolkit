# Core library
file(GLOB_RECURSE LIB_SOURCES CONFIGURE_DEPENDS
    "src/*.cpp"
    "src/*.h"
    "src/*.hpp"
)
add_library(wmtk_def_lib "${LIB_SOURCES}")
add_library(wmtk::def_lib ALIAS wmtk_def_lib)

# Include headers
target_include_directories(wmtk_def_lib PUBLIC src)

# Direct dependencies
target_link_libraries(wmtk_def_lib PUBLIC
	wmtk::toolkit
	wmtk::data
	igl::predicates
	CLI11::CLI11
#	wmtk::sec_lib
)

add_subdirectory(app)
add_subdirectory(tests)