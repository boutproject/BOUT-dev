bout_add_integrated_test(
  test-snes-save-jacobian
  SOURCES test_snes_save_jacobian.cxx
  CONFLICTS BOUT_USE_METRIC_3D
  USE_DATA_BOUT_INP
  REQUIRES BOUT_HAS_FFTW
)
