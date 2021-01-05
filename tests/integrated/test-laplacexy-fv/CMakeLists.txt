bout_add_integrated_test(test-laplacexy-fv
  SOURCES test-laplacexy.cxx
  REQUIRES BOUT_HAS_PETSC
  CONFLICTS BOUT_USE_METRIC_3D # Test requires Laplace_perpXY operator, which is not implemented for 3d metrics
  USE_RUNTEST
  USE_DATA_BOUT_INP
  )
