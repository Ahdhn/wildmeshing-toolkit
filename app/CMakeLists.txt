
include(wmtk_data)

option(WMTK_APP_TETWILD "tetwild application" ON)
option(WMTK_APP_ADAPTIVE_TESSELLATION "adaptive tessellation application" ON)
option(WMTK_APP_HARMTET "Harmonic Triangulation app" ON)
option(WMTK_APP_UNIFORM_REMESH "uniform remesh for surfaces" ON)
option(WMTK_APP_SEC "surface shortest edge collapse" ON)
option(WMTK_APP_QSLIM "surface qslim simplification" ON)
option(WMTK_APP_UNIT_TESTS "unit tests for applications" ON)
option(WMTK_APP_INTERIOR_TET_OPT "interior tet opt application" ON)
option(WMTK_APP_REPLAYER "replaying existing mesh deformation application" OFF)

if(WMTK_APP_REPLAYER)
add_subdirectory(replayer)
endif()

if(WMTK_APP_INTERIOR_TET_OPT)
add_subdirectory(interior_tet_opt)
endif()

if(WMTK_APP_TETWILD)
    add_subdirectory(tetwild)
endif()

if(WMTK_APP_ADAPTIVE_TESSELLATION)
    add_subdirectory(adaptive_tessellation)
endif()

if(WMTK_APP_HARMTET)
    add_subdirectory(harmonic_tet)
endif()

if(WMTK_APP_UNIFORM_REMESH)
    add_subdirectory(remeshing)
endif()

if(WMTK_APP_SEC)
    add_subdirectory(shortest_edge_collapse)
endif()

if(WMTK_APP_QSLIM)
    add_subdirectory(qslim)
endif()

