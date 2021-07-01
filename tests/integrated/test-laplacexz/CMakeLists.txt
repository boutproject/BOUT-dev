bout_add_integrated_test(test-laplacexz
  SOURCES test-laplacexz.cxx
  REQUIRES BOUT_HAS_PETSC
  CONFLICTS BOUT_USE_METRIC_3D
  USE_RUNTEST
  USE_DATA_BOUT_INP)
