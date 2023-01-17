bout_add_integrated_test(test-boutpp-simple-model
  USE_RUNTEST
  EXTRA_FILES mini/BOUT.inp
  REQUIRES BOUT_ENABLE_PYTHON
)

add_dependencies(test-boutpp-simple-model boutpp)
add_dependencies(build-check-integrated-tests boutpp)
