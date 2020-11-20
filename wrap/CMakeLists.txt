#set(python_module_name _brille)
#pybind11_add_module(${python_module_name} MODULE _brille.cpp "")

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
  #ensure that the python interpreter knows how to import numpy
  execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c "import numpy"
    RESULT_VARIABLE EXIT_CODE
    OUTPUT_QUIET
  )
  if (${EXIT_CODE} EQUAL 0)
    # python interpreter knows about numpy, so we can add all test_*.py tests
    set(test_folder "${CMAKE_CURRENT_SOURCE_DIR}/tests")
    file(GLOB python_tests_ tests/test_*.py)
    foreach(python_test ${python_tests_})
      add_test(NAME ${python_test} COMMAND ${PYTHON_EXECUTABLE} ${python_test} )
    endforeach()
  else()
    message(WARNING "Install numpy for ${PYTHON_EXECUTABLE} to run all tests")
  endif()
endif()
