# Directory specific compilation flags
if (MSVC)
    # warning level 4 -- add /WX for all warnings as errors
    # /Zc__cplusplus ensures the preprocessor macro __cplusplus is set
    add_compile_options(/W4 /Zc__cplusplus)
else()
    # lots of warnings -- add -Werror for  all warnings as errors
    add_compile_options(-Wall -Wextra -pedantic -Wmissing-braces)
endif()

# if(CMAKE_CXX_CPPCHECK)
#   list(APPEND CMAKE_CXX_CPPCHECK
#     "--suppress=compareBoolExpressionWithInt:${CMAKE_CURRENT_LIST_DIR}/arrayvector.hpp:349"
#   )
# endif()


set(CXX_SOURCES
  bravais.cpp
  bz.cpp
  bz_wedge.cpp
  comparisons.cpp
  debug.cpp
  hall_symbol.cpp
  lattice.cpp
  neighbours.cpp
  pointgroup.cpp
  pointsymmetry.cpp
  spg_database.cpp
  symmetry.cpp
)

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
