set(BRILLE_PYTHON_MODULE_SOURCES
  "${BRILLE_PYTHON_MODULE}.cpp" # this needs to be first
  _approx.cpp
  _basis.cpp
  _bravais.cpp
  _bz.cpp
  _debug.cpp
  _enums.cpp
  _hall_symbol.cpp
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

# Only run Python tests on the installed module

#if(BRILLE_BUILD_TESTING)
#  #ensure that the python interpreter knows how to import numpy
#  execute_process(
#    COMMAND ${Python3_EXECUTABLE} -c "import numpy"
#    RESULT_VARIABLE EXIT_CODE
#    OUTPUT_QUIET
#  )
#  if (${EXIT_CODE} EQUAL 0)
#    # python interpreter knows about numpy, so we can add all test_*.py tests
#    set(test_folder "${CMAKE_CURRENT_SOURCE_DIR}/tests")
#    file(GLOB python_tests_ tests/test_*.py)
#    foreach(python_test ${python_tests_})
#      add_test(NAME ${python_test} COMMAND ${Python3_EXECUTABLE} ${python_test} )
#    endforeach()
#  else()
#    message(WARNING "Install working numpy for ${Python3_EXECUTABLE} to run all tests")
#  endif()
#endif()
