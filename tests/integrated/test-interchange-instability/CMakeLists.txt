bout_add_integrated_test(test-interchange-instability
  SOURCES 2fluid.cxx
  CONFLICTS BOUT_USE_METRIC_3D
  USE_RUNTEST
  EXTRA_FILES slab.6b5.r1.cdl slab.6b5.r10.cdl data_1/BOUT.inp data_10/BOUT.inp
  REQUIRES BOUT_HAS_NETCDF
  )
