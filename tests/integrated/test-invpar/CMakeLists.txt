bout_add_integrated_test(test-invpar
  SOURCES test_invpar.cxx
  REQUIRES NOT BOUT_USE_METRIC_3D
  USE_RUNTEST
  USE_DATA_BOUT_INP
  )
