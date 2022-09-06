bout_add_integrated_test(test-boutcore-print
  USE_RUNTEST
  EXTRA_FILES test/BOUT.inp test.py
  REQUIRES BOUT_ENABLE_PYTHON
)
