bout_add_mms_test(MMS-time-petsc
  USE_RUNTEST
  REQUIRES
    BOUT_HAS_PETSC
    BOUT_RUN_ALL_TESTS
  EXTRA_DEPENDS MMS-time
)
