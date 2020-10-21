bout_add_integrated_test(test-multigrid-laplace
  SOURCES test_multigrid_laplace.cxx
  CONFLICTS BOUT_USE_METRIC_3D
  USE_RUNTEST
  USE_DATA_BOUT_INP
  )
