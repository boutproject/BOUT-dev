bout_add_integrated_test(
  test-multigrid-laplace
  SOURCES test_multigrid_laplace.cxx
  CONFLICTS BOUT_USE_METRIC_3D
  USE_RUNTEST
  USE_DATA_BOUT_INP
  PROCESSORS 3
)

bout_copy_file(data/BOUT_jy4.inp)
bout_copy_file(data/BOUT_jy63.inp)
bout_copy_file(data/BOUT_jy127.inp)
bout_copy_file(data/BOUT_unsheared.inp)