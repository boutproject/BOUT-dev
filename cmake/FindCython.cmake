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
  RESULT_VARIABLE CYTHON_FOUND
  OUTPUT_VARIABLE CYTHON_VERSION
  )

if (CYTHON_FOUND EQUAL 0)
  set(CYTHON_FOUND ON)
else()
  set(CYTHON_FOUND OFF)
endif()
