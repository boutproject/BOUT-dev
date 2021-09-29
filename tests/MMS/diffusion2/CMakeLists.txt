bout_add_mms_test(MMS-diffusion2
  SOURCES diffusion.cxx
  EXECUTABLE_NAME cyto
  USE_RUNTEST
  EXTRA_FILES
    X/BOUT.inp
    Y/BOUT.inp
    Z/BOUT.inp
  REQUIRES BOUT_RUN_ALL_TESTS
  )
