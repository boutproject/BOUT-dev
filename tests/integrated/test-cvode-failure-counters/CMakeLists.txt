# solver:cvode_failure_counters requires SUNDIALS >= 7. BOUT_DOWNLOAD_SUNDIALS
# builds SUNDIALS 7
if(BOUT_DOWNLOAD_SUNDIALS OR SUNDIALS_VERSION VERSION_GREATER_EQUAL 7)
  set(BOUT_CVODE_FAILURE_COUNTERS_SUPPORTED ON)
endif()

bout_add_integrated_test(
  test-cvode-failure-counters
  SOURCES test_cvode_failure_counters.cxx
  REQUIRES BOUT_HAS_CVODE BOUT_CVODE_FAILURE_COUNTERS_SUPPORTED
  USE_DATA_BOUT_INP
)
