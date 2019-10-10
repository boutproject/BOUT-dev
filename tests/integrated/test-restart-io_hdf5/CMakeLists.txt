bout_add_integrated_test(test-restart-io-hdf5
  SOURCES test-restart-io.cxx
  USE_RUNTEST
  USE_DATA_BOUT_INP
  REQUIRES BOUT_HAS_HDF5
  )
