bout_add_integrated_test(test-backtrace
  SOURCES boutexcept.cxx
  REQUIRES BOUT_USE_BACKTRACE
  CONFLICTS CMAKE_SYSTEM_NAME MATCHES "FreeBSD"
  USE_RUNTEST
  )
