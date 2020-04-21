bout_add_integrated_test(test-petsc-laplace-MAST-grid
  SOURCES test_petsc_laplace_MAST_grid.cxx
  REQUIRES BOUT_HAS_PETSC
  REQUIRES NOT DEFINED ENV{TRAVIS}
  USE_RUNTEST
  USE_DATA_BOUT_INP
  EXTRA_FILES
    grids/grid_MAST_SOL_jyis2.nc
    grids/grid_MAST_SOL_jyis34.nc
    grids/grid_MAST_SOL_jyis65.nc
    grids/grid_MAST_SOL_jyis81.nc
    grids/grid_MAST_SOL_jyis113.nc
  )
