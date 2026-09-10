bout_add_integrated_test(
  test-petsc-ordering
  SOURCES test_petsc_ordering.cxx
  REQUIRES BOUT_HAS_PETSC
  USE_RUNTEST USE_DATA_BOUT_INP
  PROCESSORS 4
)
