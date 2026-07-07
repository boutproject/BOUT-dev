bout_add_integrated_test(
  test_petsc_operators
  SOURCES test_petsc_operators.cxx
  USE_RUNTEST USE_DATA_BOUT_INP
  REQUIRES BOUT_HAS_PETSC
  REQUIRES BOUT_USE_METRIC_3D
  EXTRA_FILES rotating-ellipse.py test_petsc_operators.py
)
