bout_add_mms_test(test-fci-boundary
  SOURCES get_par_bndry.cxx
  USE_RUNTEST
  USE_DATA_BOUT_INP
  PROCESSORS 1
  )


set(gridfile ${CMAKE_CURRENT_BINARY_DIR}/grid.fci.nc)
add_custom_command(OUTPUT ${gridfile}
  COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${BOUT_PYTHONPATH}:$ENV{PYTHONPATH} ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/grid.py ${gridfile}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/../../../tools/pylib/boutconfig/__init__.py
  DEPENDS grid.py
  IMPLICIT_DEPENDS  ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Creating test-fci-boundary grid file"
  )
add_custom_target(test-fci-boundary-grid DEPENDS ${gridfile})
add_dependencies(test-fci-boundary
  test-fci-boundary-grid)
