bout_add_integrated_test(test-boutpp-mms-ddz
  USE_RUNTEST
  USE_DATA_BOUT_INP
  REQUIRES BOUT_ENABLE_PYTHON
)

add_dependencies(test-boutpp-mms-ddz boutpp)
add_dependencies(build-check-integrated-tests boutpp)
