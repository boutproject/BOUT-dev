bout_add_integrated_test(test-cyclic
  SOURCES test_cyclic.cxx
  USE_RUNTEST
  USE_DATA_BOUT_INP
  EXTRA_FILES test_io.grd.nc
  REQUIRES BOUT_HAS_NETCDF
  )
