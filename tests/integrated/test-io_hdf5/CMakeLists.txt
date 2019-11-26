bout_add_integrated_test(test-io-hdf5
  SOURCES test_io_hdf5.cxx
  USE_RUNTEST
  USE_DATA_BOUT_INP
  EXTRA_FILES test_io.grd.hdf5 data/benchmark.out.0.hdf5
  REQUIRES BOUT_HAS_HDF5
)
