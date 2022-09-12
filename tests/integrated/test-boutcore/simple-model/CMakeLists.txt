bout_add_integrated_test(test-boutcore-simple-model
  USE_RUNTEST
  EXTRA_FILES mini/BOUT.inp
  REQUIRES BOUT_ENABLE_PYTHON
)

add_dependencies(test-boutcore-simple-model boutcore)
