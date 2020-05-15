bout_add_integrated_test(test-snb
  SOURCES test_snb.cxx
  REQUIRES NOT BOUT_USE_METRIC_3D
  USE_DATA_BOUT_INP
  )
