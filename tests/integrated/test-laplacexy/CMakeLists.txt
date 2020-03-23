bout_add_integrated_test(test-laplacexy
  SOURCES test-laplacexy.cxx
  REQUIRES BOUT_HAS_PETSC
  USE_RUNTEST
  USE_DATA_BOUT_INP
  )
