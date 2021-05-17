# FindNumpy
# ----------
#
# Find Numpy
#
# This module will define the following variables:
#
# ::
#
#   NUMPY_FOUND - true if FFTW was found
#   NUMPY_VERSION - Location of the FFTW includes
#   NUMPY_HEADER - Required libraries


find_package(Python 3.6 REQUIRED COMPONENTS Interpreter Development)

execute_process(COMMAND ${Python_EXECUTABLE} -c "import numpy ; print(numpy.__version__)"
  RESULT_VARIABLE NUMPY_FOUND
  OUTPUT_VARIABLE NUMPY_VERSION
  )
execute_process(COMMAND ${Python_EXECUTABLE} -c "import numpy ; print(numpy.get_include())"
  RESULT_VARIABLE NUMPY_FOUND2
  OUTPUT_VARIABLE NUMPY_HEADER
  )

if (NUMPY_FOUND EQUAL 0 AND NUMPY_FOUND2 EQUAL 0)
  set(NUMPY_FOUND ON)
else()
  set(NUMPY_FOUND OFF)
endif()
