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

execute_process(COMMAND ${Python_EXECUTABLE} -c "import jinja2"
  RESULT_VARIABLE jinja2_FOUND)
if (jinja2_FOUND EQUAL 0)
  # We have jinja2 - all good
else()
  bout_python_maybe_error(OFF jinja2)
endif()

execute_process(COMMAND ${Python_EXECUTABLE} -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX')[:-3])"
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
set(files "boutpp.pyx" "resolve_enum.pxd" "helper.cxx" "helper.h" "boutcpp.pxd")
foreach(file IN LISTS files)
  # helper.py and resolve_enum_inv.pyx.in are only required by boutpp.pyx
  #set(deps {src}/$file.in ${src}/common.sh)
  #if (${file} STREQUAL boutpp.pyx)
  #list(APPEND deps
  set(gen ${tar}/${file})
  list(APPEND generated ${gen})
  #message(FATAL_ERROR "${gen} ${src}/${file}.jinja")
  add_custom_command(OUTPUT ${gen}
	COMMAND ${CMAKE_COMMAND} -E make_directory ${tar}
	COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${tar}/..:\${PYTHONPATH} ${Python_EXECUTABLE} generate.py ${file}.jinja ${gen}
	DEPENDS ${src}/${file}.jinja
	DEPENDS ${src}/helper.py
	DEPENDS ${src}/resolve_enum_inv.pyx.jinja
	DEPENDS ${src}/generate.py
	DEPENDS bout++
	WORKING_DIRECTORY ${src}/
	COMMENT "Generating ${file}")
endforeach()

set(boutpp_depends ${generated})

set(files "boutexception_helper.hxx" "boutexception_helper.cxx" "boutpp_openmpi_compat.hxx" "bout_options.pxd" "setup.py")
foreach(file IN LISTS files)
  list(APPEND ${boutpp_depends} "${CMAKE_CURRENT_BINARY_DIR}/${file}")
  bout_copy_file("${file}")
endforeach()

add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/libboutpp.cpp
  COMMAND ${CMAKE_COMMAND} -E copy boutpp.pyx libboutpp.pyx
  COMMAND ${Python_EXECUTABLE} -m cython libboutpp.pyx --cplus -3  -X binding=True -X embedsignature=True
  COMMENT "Cythonizing python interface"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS ${boutpp_depends}
  )

add_library(boutpp${PYTHON_EXT_SUFFIX} SHARED
  ${tar}/libboutpp.cpp
  ${tar}/helper.cxx
  ${tar}/boutexception_helper.cxx
  )

add_custom_target(boutpp ALL
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/libboutpp${PYTHON_EXT_SUFFIX}.so  ${CMAKE_CURRENT_BINARY_DIR}/../boutpp/libboutpp${PYTHON_EXT_SUFFIX}.so
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/boutpp.py ${CMAKE_CURRENT_BINARY_DIR}/../boutpp/__init__.py
  DEPENDS boutpp${PYTHON_EXT_SUFFIX}
  COMMENT "Building python interface"
)

install(TARGETS  boutpp${PYTHON_EXT_SUFFIX}
  DESTINATION ${CMAKE_INSTALL_PYTHON_SITEARCH}/boutpp/
)

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/boutpp.py
  DESTINATION ${CMAKE_INSTALL_PYTHON_SITEARCH}/boutpp/
  RENAME __init__.py
  )
target_link_libraries(boutpp${PYTHON_EXT_SUFFIX} bout++)
target_include_directories(boutpp${PYTHON_EXT_SUFFIX} PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> ${Numpy_INCLUDE_DIRS} ${Python_INCLUDE_DIRS})
