bout_add_integrated_test(test-invpar
  SOURCES test_invpar.cxx
  EXECUTABLE_NAME test_invpar
  USE_RUNTEST
  USE_DATA_BOUT_INP
  EXTRA_FILES test_io.grd.nc
  )
