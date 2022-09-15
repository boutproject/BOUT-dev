bout_add_integrated_test(test-boutcore-collect-staggered
  USE_RUNTEST
  USE_DATA_BOUT_INP
  REQUIRES BOUT_ENABLE_PYTHON
)

add_dependencies(test-boutcore-collect-staggered boutcore)
add_dependencies(build-check-integrated-tests boutcore)
