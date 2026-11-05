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

find_package(Python3 3.6 COMPONENTS Interpreter Development)

if(NOT Python3_FOUND)
  message(
    STATUS
      "Could not find numpy as python was not found. Maybe the developement package is missing?"
  )
  set(Numpy_FOUND ${Python3_FOUND})
  return()
endif()

# Query the selected interpreter on every configure: the Python environment or
# its NumPy installation may have changed since the include path was cached.
execute_process(
  COMMAND "${Python3_EXECUTABLE}" -c "import numpy ; print(numpy.__version__)"
  RESULT_VARIABLE _numpy_version_result
  OUTPUT_STRIP_TRAILING_WHITESPACE
  OUTPUT_VARIABLE Numpy_VERSION
  ERROR_QUIET
)
execute_process(
  COMMAND "${Python3_EXECUTABLE}" -c "import numpy ; print(numpy.get_include())"
  RESULT_VARIABLE _numpy_include_result
  OUTPUT_STRIP_TRAILING_WHITESPACE
  OUTPUT_VARIABLE _numpy_include_dirs
  ERROR_QUIET
)

if(Numpy_DEBUG)
  message(
    STATUS
      "Looking for numpy headers in: ${_numpy_include_dirs} ${Python3_INCLUDE_DIRS}"
  )
endif()

# Do not fall back to system headers: Cython uses this interpreter's NumPy
# declarations, which must match the headers used to compile the generated code.
if(_numpy_version_result STREQUAL "0"
   AND _numpy_include_result STREQUAL "0"
   AND EXISTS "${_numpy_include_dirs}/numpy/arrayobject.h")
  set(Numpy_INCLUDE_DIR "${_numpy_include_dirs}")
else()
  set(Numpy_INCLUDE_DIR "Numpy_INCLUDE_DIR-NOTFOUND")
  unset(Numpy_VERSION)
endif()
set(Numpy_INCLUDE_DIR "${Numpy_INCLUDE_DIR}" CACHE PATH "NumPy include directory" FORCE)

if(NOT Numpy_INCLUDE_DIR)
  message(
    STATUS
      "Numpy headers not found -- do you need to install the development package?"
  )
endif()

set(Numpy_INCLUDE_DIRS ${Numpy_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  Numpy
  VERSION_VAR Numpy_VERSION
  REQUIRED_VARS Numpy_INCLUDE_DIR
)

mark_as_advanced(Numpy_INCLUDE_DIR)
