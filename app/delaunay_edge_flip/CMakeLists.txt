# Core library (header-only)
add_library(wmtk_def_lib INTERFACE)
add_library(wmtk::def_lib ALIAS wmtk_def_lib)

# Include headers
target_include_directories(wmtk_def_lib INTERFACE src)

# Direct dependencies
target_link_libraries(wmtk_def_lib INTERFACE
		wmtk::toolkit
		wmtk::data
		igl::predicates
		CLI11::CLI11
		# wmtk::sec_lib
		)

add_subdirectory(app)
add_subdirectory(tests)