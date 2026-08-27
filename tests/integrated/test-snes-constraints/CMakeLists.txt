bout_add_integrated_test(
  test-snes-constraints
  SOURCES test_snes_constraints.cxx
  REQUIRES BOUT_HAS_PETSC
  USE_DATA_BOUT_INP
)
