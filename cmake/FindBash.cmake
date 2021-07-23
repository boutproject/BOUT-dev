# FindBash
# ----------
#
# Find Bash
#
# This module will define the following variables:
#
# ::
#
#   Bash_FOUND - true if Bash was found
#   Bash_VERSION - Bash version
#   Bash_EXECUTABLE - Path to bash executable

find_program(Bash_EXECUTABLE
  bash
  )

mark_as_advanced(Bash_EXECUTABLE)

if (Bash_EXECUTABLE)
  execute_process(COMMAND "${Bash_EXECUTABLE}" --version
    RESULT_VARIABLE _bash_runs
    OUTPUT_VARIABLE _bash_stdout
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  if(_bash_stdout MATCHES "version ([0-9]+\\.[0-9]+\\.[0-9]+)")
    set(Bash_VERSION "${CMAKE_MATCH_1}")
  else()
    message (WARNING "Failed to determine version of Bash interpreter (${Bash_EXECUTABLE})! Error:\n${_Bash_STDERR}")
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Bash
  VERSION_VAR Bash_VERSION
  REQUIRED_VARS Bash_EXECUTABLE
  )
