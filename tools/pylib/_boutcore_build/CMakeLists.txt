set(BOUT_ENABLE_PYTHON MAYBE CACHE STRING "Build the Python interface")

if (BOUT_PYTHON)
  find_package(Numpy)

  if (NOT ${NUMPY_FOUND})
    set(MSG "NumPy is required for the python interface!")
    if ("${BOUT_PYTHON}" STREQUAL "ON")
      message(FATAL_ERROR ${MSG})
    else()
      message(WARNING ${MSG})
      set(BOUT_PYTHON OFF)
    endif()
  endif()

  find_package(Cython)
  if (NOT ${CYTHON_FOUND})
    set(MSG "Cython is required for the python interface!")
    if ("${BOUT_PYTHON}" STREQUAL "ON")
      message(FATAL_ERROR ${MSG})
    else()
      message(WARNING ${MSG})
      set(BOUT_PYTHON OFF)
    endif()
  endif()

  execute_process(COMMAND bash --version
    RESULT_VARIABLE BASH_FOUND
    OUTPUT_VARIABLE BASH_VERSION_INFO
    TIMEOUT 1
    )
  if (NOT ${BASH_FOUND} EQUAL 0)
    set(MSG "bash is required for the python interface!")
    if ("${BOUT_PYTHON}" STREQUAL "ON")
      message(FATAL_ERROR ${MSG})
    else()
      message(WARNING ${MSG})
      set(BOUT_PYTHON OFF)
    endif()
  endif()

  execute_process(COMMAND ${Python_EXECUTABLE} -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))"
    RESULT_VARIABLE PYTHON_WORKING
    OUTPUT_VARIABLE PYTHON_EXT_SUFFIX
    )
  if (NOT ${PYTHON_WORKING} EQUAL 0)
    set(MSG "Failed to get the extension name from python!")
    if ("${BOUT_PYTHON}" STREQUAL "ON")
      message(FATAL_ERROR ${MSG})
    else()
      message(WARNING ${MSG})
      set(BOUT_PYTHON OFF)
    endif()
  endif()
  string(STRIP ${PYTHON_EXT_SUFFIX} PYTHON_EXT_SUFFIX)

  if (NOT BUILD_SHARED_LIBS)
    set(MSG "Need to build as shared lib for the python interface!")
    if ("${BOUT_PYTHON}" STREQUAL "ON")
      message(FATAL_ERROR ${MSG})
    else()
      message(WARNING ${MSG})
      set(BOUT_PYTHON OFF)
    endif()
  endif()

  # No errors? We can build the interface!
  if ("${BOUT_PYTHON}" STREQUAL "MAYBE")
    set(BOUT_PYTHON ON)
  endif()


  if (NOT "${BOUT_PYTHON}" STREQUAL "ON")
    message(WARNING "Python interface will not be build, see warnings above")
  else()
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
      add_custom_command(OUTPUT ${gen}
	COMMAND ${CMAKE_COMMAND} -E make_directory ${tar}
	COMMAND ${CMAKE_COMMAND} -E copy ${src}/${file} ${gen}
	DEPENDS ${src}/${file}
	)
    endforeach()

    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/boutcore.cpp
      COMMAND ${Python_EXECUTABLE} -m cython boutcore.pyx --cplus -3  -X binding=True -X embedsignature=True
      COMMENT "Cythonizing python interface"
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      DEPENDS ${generated}
      )

    ADD_LIBRARY(boutcore_ SHARED
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
    target_link_libraries(boutcore_ bout++ ${Python_LIBRARIES})
    target_include_directories(boutcore_ PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> ${NUMPY_HEADERS} ${Python_INCLUDE_DIRS})

  endif() # Build instructions python interface

  #message(FATAL_ERROR "All ${PYTHON_EXT_SUFFIX} good")
endif() # Check for python interface
