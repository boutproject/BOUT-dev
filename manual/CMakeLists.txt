# BOUT++ Documentation

find_package(Doxygen)

find_package(Sphinx REQUIRED)
set(BOUT_SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/sphinx)
set(BOUT_SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/docs)

add_custom_target(sphinx-html
  COMMAND ${SPHINX_EXECUTABLE} -b html ${BOUT_SPHINX_SOURCE} ${BOUT_SPHINX_BUILD}
  COMMAND ${CMAKE_COMMAND} -E echo "Generated HTML docs in file://${BOUT_SPHINX_BUILD}/index.html"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating HTML documentation with Sphinx in ${BOUT_SPHINX_BUILD}"
)

add_custom_target(sphinx-pdf
  COMMAND ${SPHINX_EXECUTABLE} -M latexpdf ${BOUT_SPHINX_SOURCE} ${BOUT_SPHINX_BUILD}
  COMMAND ${CMAKE_COMMAND} -E echo "Generated PDF docs in file://${BOUT_SPHINX_BUILD}"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating PDF documentation with Sphinx in ${BOUT_SPHINX_BUILD}"
)

set_target_properties(sphinx-html sphinx-pdf PROPERTIES
  ENVIRONMENT PYTHONPATH=${BOUT_PYTHONPATH}:$ENV{PYTHONPATH}
)
