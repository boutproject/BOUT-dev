bout_add_integrated_test(test-naulin-laplace
  SOURCES test_naulin_laplace.cxx
  CONFLICTS BOUT_USE_METRIC_3D
  USE_RUNTEST
  USE_DATA_BOUT_INP
  REQUIRES BOUT_HAS_FFTW
  PROCESSORS 3
  )
