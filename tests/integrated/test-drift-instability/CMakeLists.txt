bout_add_integrated_test(test-drift-instability
  SOURCES 2fluid.cxx
  REQUIRES NOT BOUT_USE_METRIC_3D
  USE_RUNTEST
  USE_DATA_BOUT_INP
  EXTRA_FILES runtest.py uedge.grd_std.cdl
  )
