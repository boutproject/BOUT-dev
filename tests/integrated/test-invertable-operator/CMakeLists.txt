bout_add_integrated_test(test-invertable-operator
  SOURCES invertable_operator.cxx
  CONFLICTS BOUT_USE_METRIC_3D
  USE_RUNTEST
  USE_DATA_BOUT_INP
  REQUIRES BOUT_HAS_PETSC
  )
