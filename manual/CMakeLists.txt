# BOUT++ Documentation

find_package(Doxygen)

find_package(Sphinx REQUIRED)
set(BOUT_SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/sphinx)
set(BOUT_SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/docs)

set(env_command
    ${CMAKE_COMMAND} -E env
    PYTHONPATH=${BOUT_PYTHONPATH}:$ENV{PYTHONPATH}
)

add_custom_target(sphinx-html
  COMMAND ${env_command}
  COMMAND ${SPHINX_EXECUTABLE} -b html ${BOUT_SPHINX_SOURCE} ${BOUT_SPHINX_BUILD}
  COMMAND ${CMAKE_COMMAND} -E echo "Generated HTML docs in file://${BOUT_SPHINX_BUILD}/index.html"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating HTML documentation with Sphinx in ${BOUT_SPHINX_BUILD}"
)

add_custom_target(sphinx-pdf
  COMMAND ${env_command}
  COMMAND ${SPHINX_EXECUTABLE} -M latexpdf ${BOUT_SPHINX_SOURCE} ${BOUT_SPHINX_BUILD}
  COMMAND ${CMAKE_COMMAND} -E echo "Generated PDF docs in file://${BOUT_SPHINX_BUILD}"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating PDF documentation with Sphinx in ${BOUT_SPHINX_BUILD}"
)

add_custom_target(docs)
add_dependencies(docs sphinx-html)

install(DIRECTORY ${BOUT_SPHINX_BUILD}/
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/doc/bout++/
  EXCLUDE_FROM_ALL
  COMPONENT docs
  PATTERN .* EXCLUDE
)
