bout_add_integrated_test(test-gyro
  SOURCES test_gyro.cxx
  USE_RUNTEST
  USE_DATA_BOUT_INP
  EXTRA_FILES cyclone_68x32.nc data/benchmark.0.nc
  REQUIRES BOUT_HAS_NETCDF
  PROCESSORS 4
  )
