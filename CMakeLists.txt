cmake_minimum_required(VERSION 3.18.2)
cmake_policy(SET CMP0076 NEW) # Ensure target_sources converts relative paths
cmake_policy(SET CMP0017 NEW) # Prefer cmake's own files for include/find_package before CMAKE_MODULE_PATH

list(INSERT CMAKE_MODULE_PATH 0 "${CMAKE_CURRENT_LIST_DIR}/cmake")

project(brille)
set(BRILLE_LIBRARY_TARGET brille)
set(BRILLE_PYTHON_MODULE _brille)
set(BRILLE_PYTHON_DESTINATION brille)
set(BRILLE_SINGLE_HEADER brille.h) # must match template file in project root
SET(BRILLE_LIB_DESTINATION lib)
SET(BRILLE_BIN_DESTINATION bin)
SET(BRILLE_INCLUDE_DESTINATION include)
set(TESTING_TARGET tester)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

# Make the availability of testing optional
option(BRILLE_BUILD_TESTING "Build tests for brille" ON)
# Allow system pybind11 to be required
option(REQUIRE_SYSTEM_PYBIND11 "Never attempt to fetch pybind11" OFF)
mark_as_advanced(REQUIRE_SYSTEM_PYBIND11)
# Allow system catch2 to be required
option(REQUIRE_SYSTEM_CATCH2 "Never attempt to fetch catch2" OFF)
mark_as_advanced(REQUIRE_SYSTEM_CATCH2)
# Allow the log level to be set by cmake
set(BRILLE_LOGLEVEL "INFO" CACHE STRING "Emit log messages to standard out (DEBUG|VERBOSE)")
# Special option for profiling runs
option(BRILLE_PROFILING "Emit profiling output to standard out" OFF)

# Define the minimum version of libraries needed
set(MINIMUM_PYBIND11_VERSION 2.10.3)
set(FETCH_PYBIND11_REPO https://github.com/pybind/pybind11)
set(MINIMUM_CATCH2_VERSION 2.13.9)
set(FETCH_CATCH2_REPO https://github.com/catchorg/Catch2)

# Optional support of HDF5-based IO
option(BRILLE_HDF5 "Add HDF5 file-based IO" ON)
option(REQUIRE_SYSTEM_HIGHFIVE "Never attempt to fetch HighFive" OFF)
if (BRILLE_HDF5)
    # Include the HDF5 IO methods -- safe since configuring will fail if HDF5/HighFive is not available
    add_definitions(-DUSE_HIGHFIVE)
endif()
mark_as_advanced(REQUIRE_SYSTEM_HIGHFIVE)
set(MINIMUM_HIGHFIVE_VERSION 2.6.2)
SET(FETCH_HIGHFIVE_REPO https://github.com/BlueBrain/HighFive)

# find_program(CMAKE_CXX_CPPCHECK NAMES cppcheck)
# if (CMAKE_CXX_CPPCHECK)
#     list(
#         APPEND CMAKE_CXX_CPPCHECK
#             "--enable=warning"
#             "--inconclusive"
#             "--force"
#             "--inline-suppr"
#             # "--template=gcc" # uncomment to get suppression error ids in brackets
#     )
# endif()

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
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif(NOT CMAKE_BUILD_TYPE)

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

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
# With GCC 10+ the interproceedural optmization only adds to compilation time without improving performance
SET(CMAKE_INTERPROCEDURAL_OPTIMIZATION FALSE)

include(brille-conan)

if (PYTHON_EXECUTABLE)
  # Ensure the provided Python interpreter is used
  set(Python3_EXECUTABLE ${PYTHON_EXECUTABLE})
endif()
# With modern CMake, this find_package forces pybind11 to use FindPython instead of its custom tools
find_package(Python3 COMPONENTS Interpreter Development)

if (NOT DEFINED PYTHON_EXECUTABLE)
  set(PYTHON_EXECUTABLE ${Python3_EXECUTABLE})
endif()

include(fetcher)

# Attempt to find catch2 to handle C++ testing
if(BRILLE_BUILD_TESTING)
  git_fetch(catch2 ${MINIMUM_CATCH2_VERSION} ${FETCH_CATCH2_REPO} ${REQUIRE_SYSTEM_CATCH2})
  list(APPEND CMAKE_MODULE_PATH "${catch2_SOURCE_DIR}/contrib")
else()
  # Since no testing is to be built, fake that we've found Catch2.
  set(Catch2_FOUND ON)
endif()
# Attempt to find pybind11 to handle CPython bindings
git_fetch(pybind11 ${MINIMUM_PYBIND11_VERSION} ${FETCH_PYBIND11_REPO} ${REQUIRE_SYSTEM_PYBIND11})

# Read the version of brille
include(checkgit)
checkGitSetup(BRILLE) # defines/creates version.hpp; and library 'git_version'
# So that we can print it to the console along with the specified build type
message(STATUS "Build brille v${BRILLE_VERSION} with type ${CMAKE_BUILD_TYPE}")

if (NOT SKBUILD)
# Create a single header by contatenating all headers in src/
add_custom_target(single_header
  COMMAND ${Python3_EXECUTABLE} acme.py ${BRILLE_SINGLE_HEADER} -o ${CMAKE_BINARY_DIR}
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  DEPENDS git_version
)
add_custom_target(library DEPENDS single_header ${BRILLE_LIBRARY_TARGET})
endif()

# We will always build the python module
list(APPEND CXX_TARGETS ${BRILLE_PYTHON_MODULE})

if(BRILLE_BUILD_TESTING)
  enable_testing() # allows registration of Python tests in wrap/
endif()

# Target for python module
pybind11_add_module(${BRILLE_PYTHON_MODULE} MODULE)
add_subdirectory(wrap)

if(BRILLE_BUILD_TESTING)
  list(APPEND CXX_TARGETS ${TESTING_TARGET}) # Include the C++ test target
  # target for Catch2 based tests
  add_executable(${TESTING_TARGET})
endif()

if(NOT SKBUILD)
# target for C++ shared library
add_library(${BRILLE_LIBRARY_TARGET} src/hdf_interface.hpp)
list(APPEND CXX_TARGETS ${BRILLE_LIBRARY_TARGET})
endif()

add_subdirectory(lib)
# add the dependencies and include directories for all CXX targets:
add_subdirectory(src)  # important

foreach(CXX_TARGET IN LISTS CXX_TARGETS)
  addGitVersion(${CXX_TARGET})
endforeach()

if(BRILLE_BUILD_TESTING)
  message(STATUS "Build testing target '${TESTING_TARGET}' and configure CTest")
  target_link_libraries(${TESTING_TARGET} PUBLIC Catch2::Catch2)
  include(CTest)
  include(Catch)
  catch_discover_tests(${TESTING_TARGET})
endif()

include(brille-openmp)
include(brille-hdf5)
include(brille-doxygen)

if (NOT SKBUILD)
install(
  TARGETS ${BRILLE_LIBRARY_TARGET}
  ARCHIVE DESTINATION ${BRILLE_LIB_DESTINATION}
  LIBRARY DESTINATION ${BRILLE_LIB_DESTINATION}
  RUNTIME DESTINATION ${BRILLE_BIN_DESTINATION}
)
install(
  FILES "${CMAKE_BINARY_DIR}/${BRILLE_SINGLE_HEADER}"
  DESTINATION ${BRILLE_INCLUDE_DESTINATION}
)
endif()

install(TARGETS ${BRILLE_PYTHON_MODULE} DESTINATION ${BRILLE_PYTHON_DESTINATION})
