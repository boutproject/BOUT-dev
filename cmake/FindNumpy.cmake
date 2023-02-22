# FindNumpy
# ----------
#
# Find Numpy
#
# This module will define the following variables:
#
# ::
#
#   Numpy_FOUND
#   Numpy_VERSION
#   Numpy_INCLUDE_DIR


find_package(Python 3.6 COMPONENTS Interpreter Development)

if (NOT Python_FOUND)
  message(STATUS "Could not find numpy as python was not found. Maybe the developement package is missing?")
  set(Numpy_FOUND ${Python_FOUND})
  return()
endif()

if (NOT Numpy_FOUND)
  execute_process(COMMAND ${Python_EXECUTABLE} -c "import numpy ; print(numpy.__version__)"
    OUTPUT_STRIP_TRAILING_WHITESPACE
    OUTPUT_VARIABLE Numpy_VERSION
    )
  execute_process(COMMAND ${Python_EXECUTABLE} -c "import numpy ; print(numpy.get_include())"
    OUTPUT_STRIP_TRAILING_WHITESPACE
    OUTPUT_VARIABLE _numpy_include_dirs
    )
endif()

if (Numpy_DEBUG)
  message(STATUS "Looking for numpy headers in: ${_numpy_include_dirs} ${PYTHON_INCLUDE_DIR}")
endif()

find_path(Numpy_INCLUDE_DIR
  numpy/arrayobject.h
  PATHS "${_numpy_include_dirs}" "${PYTHON_INCLUDE_DIR}"
  PATH_SUFFIXES numpy/core/include
  )

if (NOT Numpy_INCLUDE_DIR)
  message(STATUS "Numpy headers not found -- do you need to install the development package?")
endif()

set(Numpy_INCLUDE_DIRS ${Numpy_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Numpy
  VERSION_VAR Numpy_VERSION
  REQUIRED_VARS Numpy_INCLUDE_DIR
  )

mark_as_advanced(Numpy_INCLUDE_DIR)
