bout_add_integrated_test(test-boutpp-slicing
  USE_RUNTEST
  EXTRA_FILES test/BOUT.inp test.py
  REQUIRES BOUT_ENABLE_PYTHON
)

add_dependencies(test-boutpp-slicing boutpp)
add_dependencies(build-check-integrated-tests boutpp)
