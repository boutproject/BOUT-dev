bout_add_integrated_test(test-delp2
  SOURCES test_delp2.cxx
  USE_RUNTEST
  USE_DATA_BOUT_INP
  REQUIRES BOUT_HAS_FFTW
  PROCESSORS 4
  )
