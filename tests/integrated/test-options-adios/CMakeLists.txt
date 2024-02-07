bout_add_integrated_test(test-options-adios
  SOURCES test-options-adios.cxx
  USE_RUNTEST
  USE_DATA_BOUT_INP
  REQUIRES BOUT_HAS_ADIOS
  )
