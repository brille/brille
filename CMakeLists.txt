cmake_minimum_required(VERSION 2.8.12)
project(symbz)


set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_C_FLAGS_RELEASE "-Wall -O2")
set(CMAKE_C_FLAGS_DEBUG "-g -DSPGDEBUG -DSPGWARNING")
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)


# Version numbers
file(READ ${PROJECT_SOURCE_DIR}/src/symbz.h version_file)
string(REGEX MATCH "SYMBZ_MAJOR_VERSION ([0-9]+)" symbz_major_version ${version_file})
set(symbz_major_version ${CMAKE_MATCH_1})
string(REGEX MATCH "SYMBZ_MINOR_VERSION ([0-9]+)" symbz_minor_version ${version_file})
set(symbz_minor_version ${CMAKE_MATCH_1})
string(REGEX MATCH "SYMBZ_PATCH_VERSION ([0-9]+)" symbz_patch_version ${version_file})
set(symbz_patch_version ${CMAKE_MATCH_1})
set(serial "${symbz_major_version}.${symbz_minor_version}.${symbz_patch_version}")
set(soserial "1")

# Source code
include_directories("${PROJECT_SOURCE_DIR}/src")
set(SOURCES ${PROJECT_SOURCE_DIR}/src/arithmetic.cpp
            ${PROJECT_SOURCE_DIR}/src/bz.cpp
            ${PROJECT_SOURCE_DIR}/src/lattice.cpp
            ${PROJECT_SOURCE_DIR}/src/mathfunc.cpp
            ${PROJECT_SOURCE_DIR}/src/pointgroup.cpp
            ${PROJECT_SOURCE_DIR}/src/spg_database.cpp
            ${PROJECT_SOURCE_DIR}/src/symbz.cpp
            ${PROJECT_SOURCE_DIR}/src/symmetry.cpp)

# Shared library
add_library(symbz_shared SHARED ${SOURCES})
set_property(TARGET symbz_shared PROPERTY VERSION ${serial})
set_property(TARGET symbz_shared PROPERTY SOVERSION ${soserial})
set_property(TARGET symbz_shared PROPERTY OUTPUT_NAME symbz)
if (WIN32)
	install(TARGETS symbz_shared RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
else()
	install(TARGETS symbz_shared LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
endif()

# Static link library
add_library(symbz_static STATIC ${SOURCES})
set_property(TARGET symbz_static PROPERTY VERSION ${serial})
set_property(TARGET symbz_static PROPERTY SOVERSION ${soserial})
set_property(TARGET symbz_static PROPERTY OUTPUT_NAME symbz)
install(TARGETS symbz_static ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

# Header file
install(FILES ${PROJECT_SOURCE_DIR}/src/symbz.h DESTINATION ${CMAKE_INSTALL_PREFIX}/include)

# Python c++ bindings
add_subdirectory(pybind11)
pybind11_add_module(symbzpy EXCLUDE_FROM_ALL ${PROJECT_SOURCE_DIR}/src/_symbz.cpp ${SOURCES})

# Testing via `make check`
enable_testing()
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})
set(CMAKE_CTEST_COMMAND ctest -V)
add_executable(docheck EXCLUDE_FROM_ALL ${PROJECT_SOURCE_DIR}/test/test.cpp ${SOURCES})
target_link_libraries(docheck m)
add_test(docheck docheck)
add_dependencies(check docheck)
