# BOUT++ Documentation

find_package(Doxygen)

find_package(Sphinx REQUIRED)
set(BOUT_SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/sphinx)
set(BOUT_SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/docs)

add_custom_target(sphinx-html
  COMMAND ${SPHINX_EXECUTABLE} -b html ${BOUT_SPHINX_SOURCE} ${BOUT_SPHINX_BUILD}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating HTML documentation with Sphinx in ${BOUT_SPHINX_BUILD}"
)

add_custom_target(sphinx-pdf
  COMMAND ${SPHINX_EXECUTABLE} -b pdf ${BOUT_SPHINX_SOURCE} ${BOUT_SPHINX_BUILD}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating PDF documentation with Sphinx in ${BOUT_SPHINX_BUILD}"
)
