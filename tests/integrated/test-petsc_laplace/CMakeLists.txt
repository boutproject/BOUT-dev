bout_add_integrated_test(test-petsc-laplace
  SOURCES test_petsc_laplace.cxx
  EXECUTABLE_NAME test_petsc_laplace
  REQUIRES BOUT_HAS_PETSC
  USE_RUNTEST
  USE_DATA_BOUT_INP
  )
