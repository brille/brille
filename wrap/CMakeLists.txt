#set(python_module_name _brille)
#pybind11_add_module(${python_module_name} MODULE _brille.cpp "")

# Create a header file with information about the last time tye python module was built
# This file will be updated even if no other source files are modified but it does not
# cause relinking the module (if it is the *only* file which changes)
add_custom_target(write_version_info
  COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/version_info.py ${CMAKE_BINARY_DIR}/version_info.hpp
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
add_dependencies(${BRILLE_PYTHON_MODULE} write_version_info)
#execute_process(
#  COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/version_info.py ${CMAKE_BINARY_DIR}/version_info.hpp
#  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
#)

set(BRILLE_PYTHON_MODULE_SOURCES
  "${BRILLE_PYTHON_MODULE}.cpp" # this needs to be first
  _bravais.cpp
  _bz.cpp
  _debug.cpp
  _hall_symbol.cpp
  #_interpolation_data.cpp
  _interpolator.cpp
  _lattice.cpp
  _mesh.cpp
  _nest.cpp
  _pointgroup.cpp
  _pointsymmetry.cpp
  _polyhedron.cpp
  _primitive.cpp
  _sorting_status.cpp
  _spacegroup.cpp
  _symmetry.cpp
  _trellis.cpp
)
target_sources(${BRILLE_PYTHON_MODULE} PRIVATE ${BRILLE_PYTHON_MODULE_SOURCES})

#include_directories(${CMAKE_BINARY_DIR})
target_include_directories(${BRILLE_PYTHON_MODULE} PRIVATE ${CMAKE_BINARY_DIR})

if(BRILLE_BUILD_TESTING)
  set(test_folder "${CMAKE_CURRENT_SOURCE_DIR}/tests")
  file(GLOB python_tests_ tests/test_*.py)

  foreach(python_test ${python_tests_})
    add_test(NAME python||${BRILLE_PYTHON_MODULE}||${python_test}
      COMMAND ${PYTHON_EXECUTABLE} -m unittest discover -s ${test_folder} -p ${python_test}
    )
  endforeach()

#  add_test(NAME python COMMAND ${PYTHON_EXECUTABLE} -m unittest discover -v -s ${PROJECT_SOURCE_DIR}/tests)
#  add_dependencies(check _brille)
endif()
