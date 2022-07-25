cmake_minimum_required(VERSION 3.18.2) 
cmake_policy(SET CMP0076 NEW) # Ensure target_sources converts relative paths

# Set build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE
          "Debug"
          CACHE
          STRING
          "Choose the type of build, including: Debug, Release, RelWithDebInfo, MinSizeRel, etc."
          FORCE
          )
endif()

# scikit-build claims to support FindPython, but doesn't fully?
# https://github.com/scikit-build/scikit-build/issues/506#issuecomment-706791791
if(SKBUILD)
  message(STATUS "SciKit-Build PYTHON_EXECUTABLE ${PYTHON_EXECUTABLE}")
  message(STATUS "SciKit-Build PYTHON_INCLUDE_DIR ${PYTHON_INCLUDE_DIR}")
  message(STATUS "SciKit-Build PYTHON_LIBRARY ${PYTHON_LIBRARY}")
  set(Python_EXECUTABLE "${PYTHON_EXECUTABLE}")
  set(Python_INCLUDE_DIR "${PYTHON_INCLUDE_DIR}")
  set(Python_LIBRARY "${PYTHON_LIBRARY}")
else()
  if (PYTHON_EXECUTABLE)
    # Ensure the provided Python interpreter is used
    set(Python3_EXECUTABLE ${PYTHON_EXECUTABLE})
  endif()
endif()

# Attempt to find full development Python3
find_package(Python3 COMPONENTS Interpreter Development QUIET)
if (NOT DEFINED Python3_Development_FOUND OR NOT ${Python3_Development_FOUND})
  # Fall-back to Development.Module which is present for SciKit-Build under manylinux, but not Windows?
  find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
endif()

if (NOT DEFINED PYTHON_EXECUTABLE)
  set(PYTHON_EXECUTABLE ${Python3_EXECUTABLE})
endif()

# Read the version of brille
include(checkgit.cmake)
checkGitSetup("BRILLE") # defines/creates version.hpp; and library 'git_version'
# So that we can print it to the console along with the specified build type
message(STATUS "Build brille v${BRILLE_VERSION} with type ${CMAKE_BUILD_TYPE}")

project(brille VERSION ${BRILLE_SAFE_VERSION} LANGUAGES CXX)
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
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
# With GCC 10+ the interproceedural optmization only adds to compilation time without improving performance
SET(CMAKE_INTERPROCEDURAL_OPTIMIZATION FALSE)

if (SKBUILD)
  set(BRILLE_BUILD_TESTING OFF)
  set(BRILLE_LOGLEVEL "INFO")
endif()

if (NOT DEFINED BRILLE_BUILD_TESTING)
  # Make the availability of testing optional
  option(BRILLE_BUILD_TESTING "Build tests for brille" ON)
endif()
if (NOT DEFINED BRILLE_LOGLEVEL)
  # Allow the log level to be set by cmake
  set(BRILLE_LOGLEVEL "INFO" CACHE STRING "Emit log messages to standard out (DEBUG|VERBOSE)")
endif()
# Special option for profiling runs
option(BRILLE_PROFILING "Emit profiling output to standard out" OFF)


# Use conan to fetch/configure/build pybind11, catch2 and hdf5
if(NOT EXISTS "${CMAKE_BINARY_DIR}/conan.cmake")
  message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
  file(DOWNLOAD "https://raw.githubusercontent.com/conan-io/cmake-conan/0.18.1/conan.cmake"
       "${CMAKE_BINARY_DIR}/conan.cmake"
       EXPECTED_HASH SHA256=5cdb3042632da3efff558924eecefd580a0e786863a857ca097c3d1d43df5dcd
       TLS_VERIFY ON)
endif()
include(${CMAKE_BINARY_DIR}/conan.cmake)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_BINARY_DIR})
list(APPEND CMAKE_PREFIX_PATH ${CMAKE_CURRENT_BINARY_DIR})

conan_cmake_configure(
        REQUIRES
        pybind11/2.6.2
        catch2/2.13.9
        hdf5/1.12.0
#        highfive/2.4.1
        GENERATORS
        cmake_find_package
        OPTIONS
        hdf5:shared=True
#        highfive:with_boost=False
#        highfive:with_eigen=False
#        highfive:with_xtensor=False
#        highfive:with_opencv=False
)

conan_cmake_autodetect(conan_settings)

conan_cmake_install(PATH_OR_REFERENCE ${CMAKE_CURRENT_BINARY_DIR}
        BUILD outdated
        REMOTE conancenter
        SETTINGS ${conan_settings}
        OUTPUT_QUIET)

find_package(Catch2 REQUIRED)
find_package(pybind11 REQUIRED)
#find_package(HighFive REQUIRED)
#set(HF_LIBRARY HighFive::HighFive)

include(fetcher.cmake)
# Support of HDF5-based IO via the HighFive header-only library (which doesn't work via Conan for brille)
option(REQUIRE_SYSTEM_HIGHFIVE "Never attempt to fetch HighFive" OFF)
mark_as_advanced(REQUIRE_SYSTEM_HIGHFIVE)

set(MINIMUM_HIGHFIVE_VERSION 2.4.1)
SET(FETCH_HIGHFIVE_REPO https://github.com/BlueBrain/HighFive)
set(HIGHFIVE_USE_BOOST OFF)
set(HIGHFIVE_UNIT_TESTS OFF)
set(HIGHFIVE_EXAMPLES OFF)
set(HIGHFIVE_BUILD_DOCS OFF)
git_fetch(highfive ${MINIMUM_HIGHFIVE_VERSION} ${FETCH_HIGHFIVE_REPO} ${REQUIRE_SYSTEM_HIGHFIVE})
set(HF_LIBRARY HighFive)
add_definitions(-DUSE_HIGHFIVE)

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

if (NOT SKBUILD)
# target for C++ shared library
add_library(${BRILLE_LIBRARY_TARGET} src/hdf_interface.hpp src/is_musl.h)
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

foreach(HF_TARGET IN LISTS CXX_TARGETS)
  message(STATUS "Add HighFive library to ${HF_TARGET}")
  target_link_libraries(${HF_TARGET} PUBLIC ${HF_LIBRARY})
endforeach()


# first we can indicate the documentation build as an option (default OFF)
option(BUILD_DOC "Build documentation" OFF)
option(USE_DOXYGEN "Look for and use Doxygen to build documentation" OFF)
# check if Doxygen is installed
if (USE_DOXYGEN)
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
endif()

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

install(TARGETS ${BRILLE_PYTHON_MODULE} LIBRARY DESTINATION ${BRILLE_PYTHON_DESTINATION})
