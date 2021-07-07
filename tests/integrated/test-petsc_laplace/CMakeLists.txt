bout_add_integrated_test(test-petsc-laplace
  SOURCES test_petsc_laplace.cxx
  REQUIRES BOUT_HAS_PETSC
  CONFLICTS BOUT_USE_METRIC_3D # default preconditioner uses 'cyclic' Laplace solver which is not available with 3d metrics
  USE_RUNTEST
  USE_DATA_BOUT_INP
  )
