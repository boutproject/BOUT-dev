bout_add_integrated_test(test-datafile-facade
  SOURCES test-datafile-facade.cxx
  USE_RUNTEST
  USE_DATA_BOUT_INP
  REQUIRES BOUT_HAS_NETCDF
  PROCESSORS 4
  )
