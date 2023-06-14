bout_add_mms_test(MMS-spatial-fci
  SOURCES fci_mms.cxx
  USE_RUNTEST
  USE_DATA_BOUT_INP
  REQUIRES zoidberg_FOUND
  PROCESSORS 2
)
