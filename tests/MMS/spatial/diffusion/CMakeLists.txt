bout_add_mms_test(MMS-spatial-diffusion
  SOURCES diffusion.cxx
  USE_RUNTEST
  EXTRA_FILES X/BOUT.inp
  REQUIRES BOUT_RUN_ALL_TESTS
  PROCESSORS 2
)
