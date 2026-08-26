bout_add_integrated_test(
  test-snes-save-jacobian
  SOURCES test_snes_save_jacobian.cxx
  REQUIRES BOUT_HAS_PETSC
  USE_DATA_BOUT_INP
)

configure_file(../jacobian_tools/read_jacobian.py read_jacobian.py COPYONLY)
configure_file(../jacobian_tools/PetscBinaryIO.py PetscBinaryIO.py COPYONLY)
