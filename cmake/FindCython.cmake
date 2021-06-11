# FindCython
# ----------
#
# Find Cython
#
# This module will define the following variables:
#
# ::
#
#   CYTHON_FOUND - true if Cython was found
#   CYTHON_VERSION - Cython version

execute_process(COMMAND ${Python_EXECUTABLE} -c "import cython ; print(cython.__version__)"
  RESULT_VARIABLE _cython_runs
  OUTPUT_VARIABLE CYTHON_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )

if (${_cython_runs} EQUAL 0)
  set(CYTHON_RUNS TRUE)
else()
  set(CYTHON_RUNS FALSE)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Cython
  VERSION_VAR CYTHON_VERSION
  REQUIRED_VARS CYTHON_RUNS
  )
