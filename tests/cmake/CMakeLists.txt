add_test(
  NAME cmake-language-standards
  COMMAND ${CMAKE_COMMAND}
          -S ${CMAKE_CURRENT_SOURCE_DIR}/language-standards
          -B ${CMAKE_CURRENT_BINARY_DIR}/language-standards
          -DBOUT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
)
add_test(
  NAME cmake-language-standards-build
  COMMAND ${CMAKE_COMMAND} --build
          ${CMAKE_CURRENT_BINARY_DIR}/language-standards
)
set_tests_properties(cmake-language-standards-build PROPERTIES
                     DEPENDS cmake-language-standards)
