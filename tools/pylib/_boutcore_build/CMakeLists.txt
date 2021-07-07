# Error if Python API was explicitly requested, otherwise just a
# warning and don't build Python API
macro(bout_python_maybe_error VAR NAME)
  if (NOT ${VAR})
    set(_error_msg "${NAME} is required for the Python interface")
    if (NOT "${BOUT_ENABLE_PYTHON}" STREQUAL "MAYBE")
      message(FATAL_ERROR ${_error_msg})
    else()
      message(WARNING ${_error_msg})
      set(BOUT_ENABLE_PYTHON OFF PARENT_SCOPE)
      return()
    endif()
  endif()
endmacro()


bout_python_maybe_error(BUILD_SHARED_LIBS "BOUT++ shared library")

find_package(Numpy)
bout_python_maybe_error(${Numpy_FOUND} Numpy)

find_package(Cython)
bout_python_maybe_error(${Cython_FOUND} Cython)

find_package(Bash)
bout_python_maybe_error(${Bash_FOUND} Bash)

execute_process(COMMAND ${Python_EXECUTABLE} -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))"
  RESULT_VARIABLE PYTHON_WORKING
  OUTPUT_VARIABLE PYTHON_EXT_SUFFIX
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )
if (NOT ${PYTHON_WORKING} EQUAL 0)
  set(MSG "Failed to get the extension name from python!")
  if ("${BOUT_ENABLE_PYTHON}" STREQUAL "ON")
    message(FATAL_ERROR ${MSG})
  else()
    message(WARNING ${MSG})
    set(BOUT_ENABLE_PYTHON OFF )
  endif()
endif()

# No errors? We can build the interface!
if ("${BOUT_ENABLE_PYTHON}" STREQUAL "MAYBE")
  set(BOUT_ENABLE_PYTHON ON PARENT_SCOPE)
endif()

if (NOT BOUT_ENABLE_PYTHON)
  message(WARNING "Python interface will not be built, see warnings above")
  return()
endif()
message(STATUS "Building BOUT++ Python interface")

set(generated)
set(src ${CMAKE_CURRENT_SOURCE_DIR})
set(tar ${CMAKE_CURRENT_BINARY_DIR})
set(files "setup.py" "boutcore.pyx" "resolve_enum.pxd" "helper.cxx" "helper.h" "boutcpp.pxd")
foreach(file IN LISTS files)
  # helper.py and resolve_enum_inv.pyx.in are only required by boutcore.pyx
  #set(deps {src}/$file.in ${src}/common.sh)
  #if (${file} STREQUAL boutcore.pyx)
  #list(APPEND deps
  set(gen ${tar}/${file})
  list(APPEND generated ${gen})
  #message(FATAL_ERROR "${gen} ${src}/${file}.in")
  add_custom_command(OUTPUT ${gen}
	COMMAND ${CMAKE_COMMAND} -E make_directory ${tar}
	COMMAND bash ${src}/${file}.in --cmake > ${gen}
	DEPENDS ${src}/${file}.in
	DEPENDS ${src}/helper.py
	DEPENDS ${src}/resolve_enum_inv.pyx.in
	DEPENDS ${src}/common.sh
	DEPENDS bout++
	WORKING_DIRECTORY ${src}/
	COMMENT "Generating ${file}")
endforeach()
set(files "boutexception_helper.hxx" "boutexception_helper.cxx" "boutcore_openmpi_compat.hxx")
foreach(file IN LISTS files)
  set(gen ${CMAKE_CURRENT_BINARY_DIR}/${file})
  list(APPEND generated ${gen})
  bout_copy_file("${file}")
endforeach()

add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/boutcore.cpp
  COMMAND ${Python_EXECUTABLE} -m cython boutcore.pyx --cplus -3  -X binding=True -X embedsignature=True
  COMMENT "Cythonizing python interface"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS ${generated}
  )

add_library(boutcore_ SHARED
  ${tar}/boutcore.cpp
  ${tar}/helper.cxx
  ${tar}/boutexception_helper.cxx
  )

add_custom_target(boutcore ALL
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/libboutcore_.so ${CMAKE_CURRENT_BINARY_DIR}/../boutcore${PYTHON_EXT_SUFFIX}
  BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/../boutcore${PYTHON_EXT_SUFFIX}
  DEPENDS boutcore_
  COMMENT "Building python interface"
  )
target_link_libraries(boutcore_ bout++)
target_include_directories(boutcore_ PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> ${Numpy_INCLUDE_DIRS} ${Python_INCLUDE_DIRS})
