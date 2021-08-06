bout_add_integrated_test(test-drift-instability
  SOURCES 2fluid.cxx
  USE_RUNTEST
  USE_DATA_BOUT_INP
  EXTRA_FILES runtest.py uedge.grd_std.cdl
  REQUIRES BOUT_HAS_NETCDF
  )
