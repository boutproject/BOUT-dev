bout_add_integrated_test(test_initial
  SOURCES test_initial.cxx
  EXTRA_FILES make_test_functions_include.py
  USE_RUNTEST
  USE_DATA_BOUT_INP
  )

execute_process(COMMAND python3 make_test_functions_include.py
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  RESULT_VARIABLE TEST_INITIAL_MAKE_TEST_FUNCTIONS_RESULT)

if(NOT TEST_INITIAL_MAKE_TEST_FUNCTIONS_RESULT EQUAL "0")
  message(FATAL_ERROR "Unable to make test_functions.cxx needed for test_initial:"
    " ${TEST_INITIAL_MAKE_TEST_FUNCTIONS_RESULT}")
endif()
