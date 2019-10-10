bout_add_integrated_test(test-laplace
  SOURCES test_laplace.cxx
  EXECUTABLE_NAME test_laplace
  EXTRA_FILES test_laplace.grd.nc data/benchmark.0.nc
  USE_RUNTEST
  USE_DATA_BOUT_INP
  )
