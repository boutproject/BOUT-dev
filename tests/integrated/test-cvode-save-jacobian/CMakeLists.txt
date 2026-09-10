bout_add_integrated_test(
  test-cvode-save-jacobian
  SOURCES test_cvode_save_jacobian.cxx
  REQUIRES BOUT_HAS_PETSC BOUT_HAS_CVODE
  USE_DATA_BOUT_INP
)

configure_file(../jacobian_tools/read_jacobian.py read_jacobian.py COPYONLY)
configure_file(../jacobian_tools/PetscBinaryIO.py PetscBinaryIO.py COPYONLY)
