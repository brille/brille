cmake_minimum_required(VERSION 3.13) # tested working 3.13.0, 3.13.4, 3.14.5
cmake_policy(SET CMP0076 NEW) # Ensure target_sources converts relative paths

project(brille)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# find_program(CMAKE_CXX_CPPCHECK NAMES cppcheck)
# if (CMAKE_CXX_CPPCHECK)
#     list(
#         APPEND CMAKE_CXX_CPPCHECK
#             "--enable=warning"
#             "--inconclusive"
#             "--force"
#             "--inline-suppr"
#             # "--template=gcc" # uncomment to get suppression error ids in brackets
#             "--suppress=*:${CMAKE_SOURCE_DIR}/lib/*"
#             "--suppress=compareBoolExpressionWithInt:${CMAKE_SOURCE_DIR}/src/arrayvector.hpp:349"
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
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(Python_ADDITIONAL_VERSIONS "3")
find_package(PythonInterp)

# Python c++ bindings
find_package(pybind11 REQUIRED)
pybind11_add_module(_brille ${PROJECT_SOURCE_DIR}/src/_brille.cpp "")

add_custom_target(write_version_info
                  COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/version_info.py ${CMAKE_BINARY_DIR}/version_info.hpp
                  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
add_dependencies(_brille write_version_info)
execute_process(COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/version_info.py ${CMAKE_BINARY_DIR}/version_info.hpp WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_BINARY_DIR})
execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from version_info import version_number; version_number()"
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE BRILLE_VERSION)
message(STATUS "Building version ${BRILLE_VERSION}")



# Testing via `make check`
enable_testing()
# to see output from failed tests, run via, e.g., CTEST_OUTPUT_ON_FAILURE=1 make check
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})
set(CMAKE_CTEST_COMMAND ctest -V)
# Catch2 based tests
add_executable(catch2 EXCLUDE_FROM_ALL "")
add_dependencies(check catch2)
add_test(NAME catch2 COMMAND catch2)
# Python unittest based tests
if (${PYTHONINTERP_FOUND})
	add_test(NAME python COMMAND ${PYTHON_EXECUTABLE} -m unittest discover -v -s ${PROJECT_SOURCE_DIR}/tests)
	add_dependencies(check _brille)
endif()

set(CXX_TARGETS _brille catch2)
foreach(CXX_TARGET ${CXX_TARGETS})
  set_property(TARGET ${CXX_TARGET} PROPERTY CXX_STANDARD 14)
endforeach(CXX_TARGET)


# OpenMP support:
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
		set(OMP_TARGETS _brille catch2)
		foreach(OMP_TARGET ${OMP_TARGETS})
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
  message(STATUS "Documentation generation requires Doxygen to be available")
endif (DOXYGEN_FOUND)

add_subdirectory(lib)
add_subdirectory(src)
add_subdirectory(tests)
