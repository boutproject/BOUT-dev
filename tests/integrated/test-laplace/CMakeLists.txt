bout_add_integrated_test(test-laplace
  SOURCES test_laplace.cxx
  CONFLICTS BOUT_USE_METRIC_3D
  EXTRA_FILES test_laplace.grd.nc data/benchmark.0.nc
  USE_RUNTEST
  USE_DATA_BOUT_INP
  )
