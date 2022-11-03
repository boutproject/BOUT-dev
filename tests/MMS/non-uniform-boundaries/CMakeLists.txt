bout_add_mms_test(MMS-non-uniform-boundaries
  USE_RUNTEST
  USE_DATA_BOUT_INP
  REQUIRES BOUT_ENABLE_PYTHON
)

add_dependencies(MMS-non-uniform-boundaries boutcore)
