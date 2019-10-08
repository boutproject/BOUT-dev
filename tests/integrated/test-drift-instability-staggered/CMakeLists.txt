bout_add_integrated_test(drift-instability-staggered
  SOURCES 2fluid.cxx
  EXECUTABLE_NAME 2fluid
  USE_RUNTEST
  USE_DATA_BOUT_INP
  EXTRA_FILES runtest.py uedge.grd_std.cdl
  )
