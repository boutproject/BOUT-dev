bout_add_integrated_test(test-laplacexy2-hypre
  SOURCES test-laplacexy.cxx
  REQUIRES BOUT_HAS_HYPRE
  CONFLICTS BOUT_USE_METRIC_3D
  USE_RUNTEST
  USE_DATA_BOUT_INP
  )
