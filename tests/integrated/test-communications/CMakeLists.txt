bout_add_integrated_test(test-communications
  SOURCES test-communications.cxx
  USE_RUNTEST
  USE_DATA_BOUT_INP
  EXTRA_FILES
    data_limiter/BOUT.inp
  )
