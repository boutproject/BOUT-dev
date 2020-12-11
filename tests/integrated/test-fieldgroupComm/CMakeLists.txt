bout_add_integrated_test(test-fieldgroupComm
  SOURCES test_fieldgroupcomm.cxx
  USE_RUNTEST
  USE_DATA_BOUT_INP
  EXTRA_FILES cyclone_68x32.nc
  REQUIRES BOUT_HAS_NETCDF
  )
