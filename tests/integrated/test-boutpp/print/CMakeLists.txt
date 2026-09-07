bout_add_integrated_test(
  test-boutpp-print
  EXTRA_FILES test/BOUT.inp test.py
  REQUIRES BOUT_ENABLE_PYTHON
)

add_dependencies(test-boutpp-print boutpp)
add_dependencies(build-check-integrated-tests boutpp)
