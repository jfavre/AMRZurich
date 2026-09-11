if (PARAVIEW_USE_PYTHON)
  file(GENERATE
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/TestAMAZEReader.py
    INPUT ${CMAKE_CURRENT_SOURCE_DIR}/TestAMAZEReader.py.in)

  set(_vtk_build_TEST_FILE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  paraview_add_test_python(
    NO_RT DIRECT_DATA
    TestAMAZEReader.py)

  set_tests_properties(Python-TestAMAZEReader
    PROPERTIES ENVIRONMENT "PARAVIEW_DATA_ROOT=${CMAKE_CURRENT_SOURCE_DIR}/Testing/Data")
endif()

