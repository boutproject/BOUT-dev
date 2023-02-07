bout_add_mms_test(test-fci-mpi
  SOURCES fci_mpi.cxx
  USE_RUNTEST
  USE_DATA_BOUT_INP
  PROCESSORS 6
  DOWNLOAD https://zenodo.org/record/7614499/files/W7X-conf4-36x8x128.fci.nc?download=1
  DOWNLOAD_NAME grid.fci.nc
)
