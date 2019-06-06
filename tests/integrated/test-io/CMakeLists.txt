bout_add_integrated_test(test_io
  SOURCES test_io.cxx
  USE_RUNTEST
  USE_DATA_BOUT_INP
  EXTRA_FILES test_io.grd.nc data/benchmark.out.0.nc
  REQUIRES NETCDF
  )
