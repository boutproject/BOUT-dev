bout_add_integrated_test(test-invertable-operator
  SOURCES invertable_operator.cxx
  USE_RUNTEST
  USE_DATA_BOUT_INP
  REQUIRES BOUT_HAS_PETSC
  )
