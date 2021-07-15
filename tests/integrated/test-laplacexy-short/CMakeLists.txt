bout_add_integrated_test(test-laplacexy-short
  SOURCES test-laplacexy.cxx
  REQUIRES BOUT_HAS_PETSC
  CONFLICTS BOUT_USE_METRIC_3D # Test uses cyclic Laplace solver as a preconditioner, which is not available with 3d metrics
  USE_RUNTEST
  USE_DATA_BOUT_INP
  )
