bout_add_integrated_test(test-smooth
  SOURCES test_smooth.cxx
  USE_RUNTEST
  USE_DATA_BOUT_INP
  EXTRA_FILES test_smooth.nc data/benchmark.0.nc
  )
