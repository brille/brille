if (MSVC)
    # warning level 4 -- add /WX for all warnings as errors
    add_compile_options(/W4)
    # suppress MSVC warning C4996 about 'localtime' vs 'localtime_s'
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
    # Allow UTF-8 identifiers https://stackoverflow.com/a/47704050
    add_compile_options("$<$<C_COMPILER_ID:MSVC>:/utf-8>")
    add_compile_options("$<$<CXX_COMPILER_ID:MSVC>:/utf-8>")
else()
    # lots of warnings -- add -Werror for  all warnings as errors
    add_compile_options(-Wall -Wextra -pedantic)
endif()

set(CMAKE_MACOSX_RPATH OFF)

if (${BRILLE_LOGLEVEL} STREQUAL "VERBOSE")
    message(STATUS "Verbose logging emitted at runtime")
    add_definitions(-DVERBOSE)
else()
    if (${BRILLE_LOGLEVEL} STREQUAL "DEBUG")
        message(STATUS "Debug logging emitted at runtime")
        add_definitions(-DDEBUG)
    else()
        message(STATUS "Informational logging emitted at runtime")
    endif()
endif()
if (BRILLE_PROFILING)
    message(STATUS "Profiling output emitted at runtime")
    add_definitions(-DPROFILING)
endif (BRILLE_PROFILING)

# Directory specific compilation flags
# Setting directory-level flags here does nothing since the *target* is defined in the parent folder.
# Instead we can specify flags *per source file*
if (MSVC)
    # warning level 4 -- add /WX for all warnings as errors
    # /Zc__cplusplus ensures the preprocessor macro __cplusplus is set
    add_compile_options(/W4 /Zc__cplusplus)
    set(SOURCE_FILE_COMPILE_OPTIONS "")
else()
    # lots of warnings -- add -Werror for  all warnings as errors
    add_compile_options(-Wall -Wextra -Wmissing-braces -Wconversion -Wsign-conversion -pedantic-errors -Werror)
    set(SOURCE_FILE_COMPILE_OPTIONS "-Wall;-Wextra;-Wmissing-braces;-Wconversion;-Wsign-conversion;-pedantic-errors;-Werrors")
endif()

# if(CMAKE_CXX_CPPCHECK)
#   list(APPEND CMAKE_CXX_CPPCHECK
#     "--suppress=compareBoolExpressionWithInt:${CMAKE_CURRENT_LIST_DIR}/arrayvector.hpp:349"
#   )
# endif()

list(APPEND CXX_SOURCES
  approx_config.cpp
  bravais.cpp
  bz.cpp
  bz_move.cpp
  bz_wedge.cpp
  comparisons.cpp
  debug.cpp
  hall_symbol.cpp
  neighbours.cpp
  polyhedron_faces.cpp
  pointgroup.cpp
  pointsymmetry.cpp
  process_id.cpp
  spg_database.cpp
  symmetry.cpp
  vertex_map_set.cpp
)

if(BRILLE_HDF5)
    list(APPEND CXX_SOURCES hdf_interface.cpp)
endif()

foreach(CXX_SOURCE IN LISTS CXX_SOURCES)
#    message(STATUS "Setting compile options for ${CXX_SOURCE} to ${SOURCE_FILE_COMPILE_OPTIONS}")
    set_source_files_properties(${CXX_SOURCE} PROPERTIES COMPILE_OPTIONS "${SOURCE_FILE_COMPILE_OPTIONS}")
endforeach()

foreach(CXX_TARGET IN LISTS CXX_TARGETS)
  target_sources(${CXX_TARGET} PRIVATE ${CXX_SOURCES})
  target_include_directories(${CXX_TARGET} PUBLIC ${CMAKE_CURRENT_LIST_DIR})
endforeach()

# # This method of forcing static linking doesn't work since MacOS uses clang.
# # Furthermore this probably can't work since python modules are required(?)
# # to be dynamically linked -- though they *can* link to static libraries.
# # So the way forward might be to statically link the C++ brille library
# # and then have the Python module use that.
# if (APPLE)
#   # force static linking for the python module's OpenMP to avoid issues
#   # trying to use the module via MATLAB's python interface
#   target_link_options(_brille PUBLIC -static-libgcc -static-libstdc++)
# endif()

if(BRILLE_BUILD_TESTING)
  add_subdirectory(tests)
endif()
