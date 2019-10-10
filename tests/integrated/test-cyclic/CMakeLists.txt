bout_add_integrated_test(test-cyclic
  SOURCES test_cyclic.cxx
  EXECUTABLE_NAME test_cylic
  USE_RUNTEST
  USE_DATA_BOUT_INP
  EXTRA_FILES test_io.grd.nc
  )
