# Adapted from
# https://github.com/lefticus/cpp_starter_project/blob/master/cmake/Sanitizers.cmake
# Public domain

function(enable_sanitizers target_name)

  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
    option(ENABLE_COVERAGE "Enable coverage reporting for gcc/clang" FALSE)
    message(STATUS "Enable coverage: ${ENABLE_COVERAGE}")

    if(ENABLE_COVERAGE)
      target_compile_options(${target_name} PUBLIC --coverage -O0 -g)
      target_link_libraries(${target_name} PUBLIC --coverage)

      find_program(fastcov_FOUND fastcov)
      message(STATUS "Looking for fastcov: ${fastcov_FOUND}")
      find_program(genhtml_FOUND genhtml)
      message(STATUS "Looking for genhtml: ${fastcov_FOUND}")

      if (fastcov_FOUND AND genhtml_FOUND)
        set(COVERAGE_NAME coverage CACHE STRING "Name of coverage output file")
        set(COVERAGE_FILE "${COVERAGE_NAME}.info")
        set(COVERAGE_MSG "Open file://${PROJECT_SOURCE_DIR}/${COVERAGE_NAME}/index.html in your browser to view coverage HTML output")

        add_custom_target(code-coverage-capture
          COMMAND
            fastcov --include "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include"
            --exclude "${CMAKE_CURRENT_SOURCE_DIR}/externalpackages"
            --lcov --process-gcno
            --output "${COVERAGE_FILE}"
          COMMAND
            genhtml --output-directory "${COVERAGE_NAME}" --demangle-cpp --legend --show-details "${COVERAGE_FILE}"
          COMMAND
            "${CMAKE_COMMAND}" -E echo ${COVERAGE_MSG}
          WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
          COMMENT "Capturing coverage information"
          BYPRODUCTS
            "${COVERAGE_FILE}"
            "${COVERAGE_NAME}/index.html"
          )

        add_custom_target(code-coverage-clean
          COMMAND
          fastcov --zerocounters
          COMMENT "Cleaning coverage information"
          )
      else()
        message(STATUS "Coverage enabled, but coverage-capture not available. Please install fastcov and lcov")
      endif()

    endif()

    set(SANITIZERS "")

    option(ENABLE_SANITIZER_ADDRESS "Enable address sanitizer" FALSE)
    if(ENABLE_SANITIZER_ADDRESS)
      list(APPEND SANITIZERS "address")
    endif()

    option(ENABLE_SANITIZER_LEAK "Enable leak sanitizer" FALSE)
    if(ENABLE_SANITIZER_LEAK)
      list(APPEND SANITIZERS "leak")
    endif()

    option(ENABLE_SANITIZER_UNDEFINED_BEHAVIOR "Enable undefined behavior sanitizer" FALSE)
    if(ENABLE_SANITIZER_UNDEFINED_BEHAVIOR)
      list(APPEND SANITIZERS "undefined")
    endif()

    option(ENABLE_SANITIZER_THREAD "Enable thread sanitizer" FALSE)
    if(ENABLE_SANITIZER_THREAD)
      if("address" IN_LIST SANITIZERS OR "leak" IN_LIST SANITIZERS)
        message(WARNING "Thread sanitizer does not work with Address and Leak sanitizer enabled")
      else()
        list(APPEND SANITIZERS "thread")
      endif()
    endif()

    option(ENABLE_SANITIZER_MEMORY "Enable memory sanitizer" FALSE)
    if(ENABLE_SANITIZER_MEMORY AND CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
      if("address" IN_LIST SANITIZERS
         OR "thread" IN_LIST SANITIZERS
         OR "leak" IN_LIST SANITIZERS)
        message(WARNING "Memory sanitizer does not work with Address, Thread and Leak sanitizer enabled")
      else()
        list(APPEND SANITIZERS "memory")
      endif()
    endif()

    list(
      JOIN
      SANITIZERS
      ","
      LIST_OF_SANITIZERS)

  endif()

  # Default value gets overridden below
  set(BOUT_USE_SANITIZERS "None" PARENT_SCOPE)

  if(LIST_OF_SANITIZERS)
    if(NOT
       "${LIST_OF_SANITIZERS}"
       STREQUAL
       "")
      set(BOUT_USE_SANITIZERS ${LIST_OF_SANITIZERS} PARENT_SCOPE)
      target_compile_options(${target_name} PUBLIC -fsanitize=${LIST_OF_SANITIZERS} -fno-omit-frame-pointer)
      target_link_options(${target_name} PUBLIC -fsanitize=${LIST_OF_SANITIZERS})
    endif()
  endif()

endfunction()
