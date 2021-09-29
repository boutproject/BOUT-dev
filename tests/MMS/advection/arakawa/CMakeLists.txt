bout_add_mms_test(MMS-advection-arakawa
  USE_RUNTEST
  USE_DATA_BOUT_INP
  EXTRA_FILES parent.py
  EXTRA_DEPENDS MMS-advection
  REQUIRES BOUT_RUN_ALL_TESTS
)
