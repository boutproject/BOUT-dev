bout_add_integrated_test(test-gyro
  SOURCES test_gyro.cxx
  REQUIRES NOT BOUT_USE_METRIC_3D
  USE_RUNTEST
  USE_DATA_BOUT_INP
  EXTRA_FILES cyclone_68x32.nc data/benchmark.0.nc
  )
