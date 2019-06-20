bout_add_integrated_test(test-slepc-solver
  SOURCES test-slepc-solver.cxx
  USE_RUNTEST
  USE_DATA_BOUT_INP
  REQUIRES BOUT_HAS_SLEPC
  )
