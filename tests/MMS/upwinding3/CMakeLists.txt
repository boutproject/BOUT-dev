bout_add_mms_test(MMS-upwinding3
  USE_RUNTEST
  USE_DATA_BOUT_INP
  REQUIRES BOUT_ENABLE_PYTHON
)

add_dependencies(MMS-upwinding3 boutpp)
