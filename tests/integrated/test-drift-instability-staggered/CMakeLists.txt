bout_add_integrated_test(test-drift-instability-staggered
  SOURCES 2fluid.cxx
  CONFLICTS BOUT_USE_METRIC_3D
  USE_RUNTEST
  USE_DATA_BOUT_INP
  EXTRA_FILES uedge.grd_std.cdl
  REQUIRES BOUT_HAS_NETCDF
  REQUIRES BOUT_HAS_FFTW
  PROCESSORS 2
  )
