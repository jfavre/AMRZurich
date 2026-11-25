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

  # Full-file regression test; SKIPs gracefully when the data file is absent
  # (override the location with the AMAZE_FULL_TEST_FILE environment variable).
  file(GENERATE
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/TestAMAZEFullFile.py
    INPUT ${CMAKE_CURRENT_SOURCE_DIR}/TestAMAZEFullFile.py.in)

  paraview_add_test_python(
    NO_RT DIRECT_DATA
    TestAMAZEFullFile.py)

  if(DEFINED ENV{AMAZE_FULL_TEST_FILE})
    set_tests_properties(Python-TestAMAZEFullFile
      PROPERTIES
        ENVIRONMENT "AMAZE_FULL_TEST_FILE=$ENV{AMAZE_FULL_TEST_FILE}"
        SKIP_REGULAR_EXPRESSION "SKIP")
  else()
    set_tests_properties(Python-TestAMAZEFullFile
      PROPERTIES SKIP_REGULAR_EXPRESSION "SKIP")
  endif()
endif()

