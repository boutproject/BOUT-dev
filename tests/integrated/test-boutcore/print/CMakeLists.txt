bout_add_integrated_test(test-boutcore-print
  USE_RUNTEST
  EXTRA_FILES test/BOUT.inp test.py
  REQUIRES BOUT_ENABLE_PYTHON
)

add_dependencies(test-boutcore-print boutcore)
add_dependencies(build-check-integrated-tests boutcore)
