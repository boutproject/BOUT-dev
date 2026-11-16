bout_add_integrated_test(
  test-snes-jacobian-reuse
  SOURCES test_snes_jacobian_reuse.cxx
  REQUIRES BOUT_HAS_PETSC
  USE_DATA_BOUT_INP
)
