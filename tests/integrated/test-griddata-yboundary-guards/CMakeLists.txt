bout_add_integrated_test(test-griddata-yboundary-guards
  SOURCES test_griddata.cxx
  USE_RUNTEST
  EXTRA_FILES
    data-doublenull-0/BOUT.inp
    data-doublenull-1/BOUT.inp
    data-doublenull-2/BOUT.inp
    data-singlenull-0/BOUT.inp
    data-singlenull-1/BOUT.inp
    data-singlenull-2/BOUT.inp
  REQUIRES BOUT_HAS_NETCDF
  )
