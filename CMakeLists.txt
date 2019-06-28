cmake_minimum_required(VERSION 2.8.12)
project(symbz)

set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_C_FLAGS_RELEASE "-Wall -O2")
if (MSVC)
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} /Zi")
else()
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -g")
endif()
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(Python_ADDITIONAL_VERSIONS "3")
find_package(PythonInterp)

if(CMAKE_BUILD_TYPE STREQUAL "Release")
	add_custom_target(write_version_info
	                  COMMAND python ${CMAKE_SOURCE_DIR}/write_version_info.py ${CMAKE_BINARY_DIR}/version_info.h
	                  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
	)
endif()
execute_process(COMMAND python ${CMAKE_SOURCE_DIR}/write_version_info.py ${CMAKE_BINARY_DIR}/version_info.h WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_BINARY_DIR})
execute_process(COMMAND python -c "from write_version_info import print_version_number; print_version_number()"
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE SYMBZ_VERSION)
message("Building for version ${SYMBZ_VERSION}")

# Source code
include_directories("${PROJECT_SOURCE_DIR}/src")
set(SOURCES ${PROJECT_SOURCE_DIR}/src/arithmetic.cpp
            ${PROJECT_SOURCE_DIR}/src/bz.cpp
            ${PROJECT_SOURCE_DIR}/src/lattice.cpp
            ${PROJECT_SOURCE_DIR}/src/neighbours.cpp
            ${PROJECT_SOURCE_DIR}/src/pointgroup.cpp
            ${PROJECT_SOURCE_DIR}/src/spg_database.cpp
            ${PROJECT_SOURCE_DIR}/src/symmetry.cpp)

# Python c++ bindings
add_subdirectory(lib/pybind11)
pybind11_add_module(_symbz ${PROJECT_SOURCE_DIR}/src/_symbz.cpp ${SOURCES})

if(MSVC)
		if(CMAKE_BUILD_TYPE STREQUAL "Debug")
			target_link_options(_symbz PUBLIC /DEBUG)
		endif()
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Release")
	add_dependencies(_symbz write_version_info)
endif()

# Testing via `make check`
enable_testing()
# to see output from failed tests, run via, e.g.,
# CTEST_OUTPUT_ON_FAILURE=1 make check
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})
set(CMAKE_CTEST_COMMAND ctest -V)
add_executable(test_cpp EXCLUDE_FROM_ALL ${PROJECT_SOURCE_DIR}/tests/test.cpp ${SOURCES})
# target_link_libraries(test_cpp symbz_shared)
add_test(NAME C++ COMMAND test_cpp)
add_dependencies(check test_cpp)

if (${PYTHONINTERP_FOUND})
	add_test(NAME python COMMAND ${PYTHON_EXECUTABLE} -m unittest discover -v -s ${PROJECT_SOURCE_DIR})
	add_dependencies(check _symbz)
endif()

file(GLOB TEST_CATCH_SOURCES ${PROJECT_SOURCE_DIR}/tests/*_test.cpp)
add_executable(test_catch EXCLUDE_FROM_ALL ${TEST_CATCH_SOURCES} ${SOURCES})
# add_executable(test_check EXCLUDE_FROM_ALL ${PROJECT_SOURCE_DIR}/tests/linalg_test.cpp ${PROJECT_SOURCE_DIR}/tests/catch_test.cpp ${SOURCES})
target_include_directories(test_catch PUBLIC "${PROJECT_SOURCE_DIR}/lib") # where catch.hpp is located
add_dependencies(check test_catch)
add_test(NAME catch2 COMMAND test_catch)


# OpenMP support:
if( CMAKE_MAJOR_VERSION GREATER 2 AND CMAKE_MINOR_VERSION GREATER 8) # new style available since 3.9
	find_package(OpenMP)
else()
	if(MSVC)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
	else()
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -pthread")
	endif()
	if(NOT TARGET OpenMP::OpenMP_CXX)
	    find_package(Threads REQUIRED)
	    add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
	    set_property(TARGET OpenMP::OpenMP_CXX
	                 PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
	    # Only works if the same flag is passed to the linker; use CMake 3.9+ otherwise (Intel, AppleClang)
	    set_property(TARGET OpenMP::OpenMP_CXX
	                 PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_CXX_FLAGS} Threads::Threads)
	endif()
endif()
if(OpenMP_CXX_FOUND)
		target_link_libraries(_symbz       PUBLIC OpenMP::OpenMP_CXX)
		target_link_libraries(test_catch   PUBLIC OpenMP::OpenMP_CXX)
endif()

if(MSVC)
	if(CMAKE_BUILD_TYPE STREQUAL "Matlab")
		find_package(Matlab)
		target_link_libraries(_symbz PUBLIC "${Matlab_ROOT_DIR}/bin/win64/libiomp5md.lib")
		set_property(TARGET _symbz APPEND PROPERTY LINK_FLAGS /nodefaultlib:vcomp)
	else()
	endif()
else()
	target_link_libraries(_symbz PUBLIC "-fopenmp")
endif()



# first we can indicate the documentation build as an option and set it to ON by default
option(BUILD_DOC "Build documentation" OFF)

# check if Doxygen is installed
find_package(Doxygen)
if (DOXYGEN_FOUND)
    # set input and output files
    set(DOXYGEN_IN ${PROJECT_SOURCE_DIR}/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
		if (BUILD_DOC)
			add_custom_target( docs ALL
	        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
	        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	        COMMENT "Generating documentation with Doxygen"
	        VERBATIM )
		else (BUILD_DOC)
	    add_custom_target( docs
	        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
	        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	        COMMENT "Generating documentation with Doxygen"
	        VERBATIM )
		endif (BUILD_DOC)
else (DOXYGEN_FOUND)
  message("Documentation generation requires Doxygen to be available")
endif (DOXYGEN_FOUND)
