bout_add_integrated_test(test-fieldgroupComm
  SOURCES test_fieldgroupcomm.cxx
  EXECUTABLE_NAME test
  USE_RUNTEST
  USE_DATA_BOUT_INP
  EXTRA_FILES cyclone_68x32.nc
  )
