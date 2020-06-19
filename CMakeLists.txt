cmake_minimum_required(VERSION 3.13) # tested working 3.13.0, 3.13.4, 3.14.5
cmake_policy(SET CMP0076 NEW) # Ensure target_sources converts relative paths

project(brille)
set(BRILLE_PYTHON_MODULE _brille)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Make the availability of testing optional
option(BRILLE_BUILD_TESTING "Build tests for brille" ON)
# Allow system pybind11 to be required
option(REQUIRE_SYSTEM_PYBIND11 "Never attempt to fetch pybind11" OFF)
mark_as_advanced(REQUIRE_SYSTEM_PYBIND11)
# Allow system catch2 to be required
option(REQUIRE_SYSTEM_CATCH2 "Never attempt to fetch catch2" OFF)
mark_as_advanced(REQUIRE_SYSTEM_CATCH2)

# Define the minimum version of libraries needed
set(MINIMUM_PYBIND11_VERSION 2.4.3)
set(FETCH_PYBIND11_REPO https://github.com/pybind/pybind11)
set(MINIMUM_CATCH2_VERSION 2.11.3)
set(FETCH_CATCH2_REPO https://github.com/catchorg/Catch2)

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
else()
    # lots of warnings -- add -Werror for  all warnings as errors
    add_compile_options(-Wall -Wextra -pedantic)
endif()

set(CMAKE_MACOSX_RPATH 1)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif(NOT CMAKE_BUILD_TYPE)


set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Use CMake v3.12+ FindPython3 to find the system-compatible interpreter *and* libraries at once
# This solves a potential issue on convoluted Windows build systems which have multiple
# python interpreters, some of which do not have build-system compatibility (I'm looking at you, msys2 python)
# In such a case FindPythonInterp might point to the wrong python interpreter and then further configuration
# steps will fail to find any (or any compatible) libraries which match the interpreter.
find_package(Python3 COMPONENTS Interpreter Development)
# Since FindPython3 and FindPythonInterp set different variables, fake the FindPythonInterp result
# to prevent it running if called by, e.g., FindPybind11
if(Python3_FOUND)
  set(PYTHONINTERP_FOUND ON)
  set(PYTHON_EXECUTABLE ${Python3_EXECUTABLE})
  set(PYTHON_VERSION_STRING "${Python3_VERSION}")
  set(PYTHON_VERSION_MAJOR ${Python3_VERSION_MAJOR})
  set(PYTHON_VERSION_MINOR ${Python3_VERSION_MINOR})
  set(PYTHON_VERSION_PATCH ${Python3_VERSION_PATCH})
endif()

# Attempt to find catch2 to handle C++ testing
if(BRILLE_BUILD_TESTING)
  if(REQUIRE_SYSTEM_CATCH2)
    find_package(Catch2 ${MINIMUM_CATCH2_VERSION} REQUIRED)
  else()
    find_package(Catch2 ${MINIMUM_CATCH2_VERSION} QUIET)
  endif()
else()
  # Since no testing is to be built, fake that we've found Catch2.
  set(Catch2_FOUND ON)
endif()
# Attempt to find pybind11 to handle CPython bindings
if(REQUIRE_SYSTEM_PYBIND11)
  find_package(pybind11 ${MINIMUM_PYBIND11_VERSION} REQUIRED)
else()
  find_package(pybind11 ${MINIMUM_PYBIND11_VERSION} QUIET)
endif()

# Load the FetchContent module only if necessary
if(NOT (Catch2_FOUND AND pybind11_FOUND))
  include(FetchContent)
endif()

# Fetch Catch2 if necessary
if(NOT Catch2_FOUND)
  message(STATUS "Fetching catch2 v${MINIMUM_CATCH2_VERSION} from ${FETCH_CATCH2_REPO}")
  FetchContent_Declare(catch2 GIT_REPOSITORY ${FETCH_CATCH2_REPO} GIT_TAG v${MINIMUM_CATCH2_VERSION})
  FetchContent_GetProperties(catch2)
  if(NOT catch2_POPULATED)
    FetchContent_Populate(Catch2)
    add_subdirectory(${catch2_SOURCE_DIR} ${catch2_BINARY_DIR})
    list(APPEND CMAKE_MODULE_PATH "${catch2_SOURCE_DIR}/contrib")
  endif()
  set(Catch2_FOUND ON)
endif()
# Fetch pybind11 if necessary
if(NOT pybind11_FOUND)
  message(STATUS "Fetching pybind11 v${MINIMUM_PYBIND11_VERSION} from ${FETCH_PYBIND11_REPO}")
  FetchContent_Declare(pybind11 GIT_REPOSITORY ${FETCH_PYBIND11_REPO} GIT_TAG v${MINIMUM_PYBIND11_VERSION})
  FetchContent_GetProperties(pybind11)
  if(NOT pybind11_POPULATED)
    FetchContent_Populate(pybind11)
    add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
  endif()
  set(pybind11_FOUND ON)
endif()

# Read the version of brille
execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from version_info import version_number; version_number()"
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE BRILLE_VERSION)
# So that we can print it to the console along with the specified build type
message(STATUS "Build brille v${BRILLE_VERSION} with type ${CMAKE_BUILD_TYPE}")


# We will always build the python module
list(APPEND CXX_TARGETS ${BRILLE_PYTHON_MODULE})

if(BRILLE_BUILD_TESTING)
  list(APPEND CXX_TARGETS test_brille) # Include the C++ test target
  enable_testing() # allows registration of Python tests in wrap/
endif()

# Target for python module
pybind11_add_module(${BRILLE_PYTHON_MODULE} MODULE)
add_subdirectory(wrap)

if(BRILLE_BUILD_TESTING)
  # target for Catch2 based tests
  add_executable(test_brille)
  target_link_libraries(test_brille PUBLIC Catch2::Catch2)
  include(CTest)
  include(Catch)
  catch_discover_tests(test_brille)
endif()
add_subdirectory(lib)

# (eventually) target for C++ shared library
# [goes here]
add_subdirectory(src)  # important

# OpenMP support:
find_package(OpenMP REQUIRED) # Change code to support missing OpenMP?
if(OpenMP_CXX_FOUND)
  foreach(OMP_TARGET IN LISTS CXX_TARGETS)
    target_link_libraries(${OMP_TARGET} PUBLIC OpenMP::OpenMP_CXX)
  endforeach(OMP_TARGET)
  if (MSVC AND MSVC_VERSION GREATER 1919)
    add_compile_options(/openmp:experimental) # this doesn't work
  endif()
endif()

# if(MSVC)
# 	if(CMAKE_BUILD_TYPE STREQUAL "Matlab")
# 		find_package(Matlab)
# 		target_link_libraries(_brille PUBLIC "${Matlab_ROOT_DIR}/bin/win64/libiomp5md.lib")
# 		set_property(TARGET _brille APPEND PROPERTY LINK_FLAGS /nodefaultlib:vcomp)
# 	endif()
# else()
# endif()

# first we can indicate the documentation build as an option (default OFF)
option(BUILD_DOC "Build documentation" OFF)
# check if Doxygen is installed
find_package(Doxygen QUIET)
if (DOXYGEN_FOUND)
  # set input and output files
  set(DOXYGEN_IN ${PROJECT_SOURCE_DIR}/Doxyfile.in)
  set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
  # request to configure the file
  configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
  if(BUILD_DOC)
    add_custom_target( docs ALL
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating documentation with Doxygen"
    VERBATIM )
  else()
    add_custom_target( docs
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generate documentation using target 'docs'"
    VERBATIM )
  endif()
else (DOXYGEN_FOUND)
  message(STATUS "Install Doxygen to build documentation")
endif (DOXYGEN_FOUND)
