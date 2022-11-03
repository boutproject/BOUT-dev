bout_add_integrated_test(test-boutpp-collect-staggered
  USE_RUNTEST
  USE_DATA_BOUT_INP
  REQUIRES BOUT_ENABLE_PYTHON
)

add_dependencies(test-boutpp-collect-staggered boutpp)
add_dependencies(build-check-integrated-tests boutpp)
