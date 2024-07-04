set(POLYSTAR_PYTHON_MODULE_SOURCES
  "${POLYSTAR_PYTHON_MODULE}.cpp" # this needs to be first
  _approx.cpp
  _debug.cpp
  _polygon.cpp
  _polyhedron.cpp
  _bitmap.cpp
  _network.cpp
  _svg.cpp
  _wire.cpp
)
target_sources(${POLYSTAR_PYTHON_MODULE} PRIVATE ${POLYSTAR_PYTHON_MODULE_SOURCES})

# Only run Python tests on the installed module

#if(POLYSTAR_BUILD_TESTING)
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
