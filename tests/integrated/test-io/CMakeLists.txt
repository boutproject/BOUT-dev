bout_add_integrated_test(test-io
  SOURCES test_io.cxx
  EXECUTABLE_NAME test_io
  USE_RUNTEST
  USE_DATA_BOUT_INP
  EXTRA_FILES test_io.grd.nc data/benchmark.out.0.nc
  REQUIRES BOUT_HAS_NETCDF
  )
