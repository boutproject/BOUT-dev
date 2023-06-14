bout_add_integrated_test(test-boutpp-collect
  USE_RUNTEST
  EXTRA_FILES input/BOUT.inp
  REQUIRES BOUT_ENABLE_PYTHON
)

add_dependencies(test-boutpp-collect boutpp)
add_dependencies(build-check-integrated-tests boutpp)
