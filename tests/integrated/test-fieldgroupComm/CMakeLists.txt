bout_add_integrated_test(test-fieldgroupComm
  SOURCES test_fieldgroupcomm.cxx
  CONFLICTS BOUT_USE_METRIC_3D
  USE_RUNTEST
  USE_DATA_BOUT_INP
  EXTRA_FILES cyclone_68x32.nc
  REQUIRES BOUT_HAS_NETCDF
  )
