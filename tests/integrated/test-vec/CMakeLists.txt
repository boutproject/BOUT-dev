bout_add_integrated_test(test-vec
  SOURCES testVec.cxx
  USE_RUNTEST
  USE_DATA_BOUT_INP
  REQUIRES BOUT_HAS_FFTW
  PROCESSORS 4
  )
