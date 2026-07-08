bout_add_integrated_test(
  test-multigrid-laplace
  SOURCES test_multigrid_laplace.cxx
  CONFLICTS BOUT_USE_METRIC_3D
  USE_RUNTEST USE_DATA_BOUT_INP
  EXTRA_FILES data/BOUT_jy4.inp data/BOUT_jy63.inp data/BOUT_jy127.inp
              data/BOUT_unsheared.inp
  PROCESSORS 3
)
