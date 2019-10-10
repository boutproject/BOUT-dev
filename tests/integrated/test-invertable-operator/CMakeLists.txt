bout_add_integrated_test(test-invertable-operator
  SOURCES invertable_operator.cxx
  EXECUTABLE_NAME invertable_operator
  USE_RUNTEST
  USE_DATA_BOUT_INP
  REQUIRES BOUT_HAS_PETSC
  )
