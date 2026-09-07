function(escape_cxx_string input output)
  set(value "${input}")
  string(REPLACE "\\" "\\\\" value "${value}")
  string(REPLACE "\"" "\\\"" value "${value}")
  string(REPLACE "\r" "" value "${value}")
  string(REPLACE "\n" "\\n\"\n\"" value "${value}")
  set(${output}
      "${value}"
      PARENT_SCOPE
  )
endfunction()

set(git_metadata_available false)
set(git_dirty false)
set(git_diff "")

if(NOT DEFINED OUTPUT_FILE)
  message(FATAL_ERROR "OUTPUT_FILE must be set")
endif()
if(NOT DEFINED SOURCE_DIR)
  message(FATAL_ERROR "SOURCE_DIR must be set")
endif()

if(NOT DEFINED GIT_EXECUTABLE OR GIT_EXECUTABLE STREQUAL "")
  find_program(GIT_EXECUTABLE git)
endif()

if(GIT_EXECUTABLE)
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" rev-parse --is-inside-work-tree
    WORKING_DIRECTORY "${SOURCE_DIR}"
    RESULT_VARIABLE git_repo_result
    OUTPUT_VARIABLE git_repo_output
    ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  if(git_repo_result EQUAL 0 AND git_repo_output STREQUAL "true")
    set(git_metadata_available true)

    execute_process(
      COMMAND "${GIT_EXECUTABLE}" rev-parse --verify HEAD
      WORKING_DIRECTORY "${SOURCE_DIR}"
      RESULT_VARIABLE git_head_result
      OUTPUT_QUIET ERROR_QUIET
    )

    if(git_head_result EQUAL 0)
      execute_process(
        COMMAND "${GIT_EXECUTABLE}" diff --no-ext-diff --binary HEAD --
        WORKING_DIRECTORY "${SOURCE_DIR}"
        RESULT_VARIABLE git_diff_result
        OUTPUT_VARIABLE git_diff
        ERROR_QUIET
      )
      if(NOT git_diff_result EQUAL 0)
        set(git_diff "")
      endif()
    else()
      execute_process(
        COMMAND "${GIT_EXECUTABLE}" diff --no-ext-diff --binary --cached --
        WORKING_DIRECTORY "${SOURCE_DIR}"
        RESULT_VARIABLE git_cached_diff_result
        OUTPUT_VARIABLE git_cached_diff
        ERROR_QUIET
      )
      if(NOT git_cached_diff_result EQUAL 0)
        set(git_cached_diff "git diff unavailable")
      endif()

      execute_process(
        COMMAND "${GIT_EXECUTABLE}" diff --no-ext-diff --binary --
        WORKING_DIRECTORY "${SOURCE_DIR}"
        RESULT_VARIABLE git_worktree_diff_result
        OUTPUT_VARIABLE git_worktree_diff
        ERROR_QUIET
      )
      if(NOT git_worktree_diff_result EQUAL 0)
        set(git_worktree_diff "")
      endif()

      set(git_diff "${git_cached_diff}")
      if(NOT git_cached_diff STREQUAL "" AND NOT git_worktree_diff STREQUAL "")
        string(APPEND git_diff "\n")
      endif()
      string(APPEND git_diff "${git_worktree_diff}")
    endif()

    if(NOT git_diff STREQUAL "")
      set(git_dirty true)
    endif()
  endif()
endif()

escape_cxx_string("${git_diff}" git_diff_escaped)

set(output_tmp "${OUTPUT_FILE}.tmp")

configure_file(
  "${CMAKE_CURRENT_LIST_DIR}/git_metadata.cxx.in" "${output_tmp}" @ONLY
)

execute_process(
  COMMAND "${CMAKE_COMMAND}" -E copy_if_different "${output_tmp}"
          "${OUTPUT_FILE}"
)
execute_process(COMMAND "${CMAKE_COMMAND}" -E rm -f "${output_tmp}")
