bout_add_integrated_test(test-squash
  SOURCES squash.cxx
  USE_RUNTEST
  USE_DATA_BOUT_INP
  REQUIRES BOUT_HAS_NETCDF
  TESTARGS "${CMAKE_CURRENT_LIST_DIR}/../../../bin"
  )
