# FindScoreP
# ----------
#
# Find the Score-P profiling tools
#
# This module will define the following variables:
#
# ::
#
#   ScoreP_FOUND - true if ScoreP was found
#   ScoreP_EXECUTABLE - Path to the ``scorep`` compiler wrapper
#   ScoreP_INCLUDE_DIRS - Location of the ScoreP includes
#
# This module will also export the ``ScoreP::ScoreP`` target.
#
# You can also set the following variables:
#
# ``ScoreP_ROOT``
#   Specify the path to the ScoreP installation to use
#
# ``ScoreP_FLAGS``
#   The flags to pass to the ``scorep`` executable as a
#   semicolon-separated list. This defaults to ``--user;--nocompiler``
#
# ``ScoreP_COMPILER_LAUNCHER``
#   The full path to the compiler wrapper plus flags as a
#   semicolon-separated list. This defaults to
#   ``/path/to/scorep;${ScoreP_FLAGS}``
#
# ``ScoreP_DEBUG``
#   Set to TRUE to get extra debugging output

find_program(ScoreP_EXECUTABLE scorep)
mark_as_advanced(ScoreP_EXECUTABLE)

get_filename_component(ScoreP_TMP "${ScoreP_EXECUTABLE}" DIRECTORY)
get_filename_component(ScoreP_EXEC_LOCATION "${ScoreP_TMP}" DIRECTORY)

if (ScoreP_DEBUG)
  message(STATUS "[ ${CMAKE_CURRENT_LIST_FILE}:${CMAKE_CURRENT_LIST_LINE} ] "
    " ScoreP_EXECUTABLE = ${ScoreP_EXECUTABLE}"
    " ScoreP_EXEC_LOCATION = ${ScoreP_EXEC_LOCATION}")
endif()


find_path(ScoreP_INCLUDE_DIRS
  NAMES scorep/SCOREP_User.h
  DOC "Score-P include directory"
  HINTS "${ScoreP_EXEC_LOCATION}"
  PATH_SUFFIXES "include"
  )
mark_as_advanced(ScoreP_INCLUDE_DIRS)

if (ScoreP_DEBUG)
  message(STATUS "[ ${CMAKE_CURRENT_LIST_FILE}:${CMAKE_CURRENT_LIST_LINE} ] "
    " ScoreP_INCLUDE_DIRS = ${ScoreP_INCLUDE_DIRS}")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ScoreP
  REQUIRED_VARS ScoreP_INCLUDE_DIRS ScoreP_EXECUTABLE
  FAIL_MESSAGE
  "Failed to find ScoreP. Please supply the path to the scorep executable using
 -DCMAKE_PROGRAM_PATH=/path/to/scorep/bin/directory"
  )

if (ScoreP_FOUND AND NOT TARGET ScoreP::ScoreP)
  set(ScoreP_FLAGS "--user;--nocompiler" CACHE STRING "ScoreP wrapper flags")
  set(ScoreP_COMPILER_LAUNCHER "${ScoreP_EXECUTABLE};${ScoreP_FLAGS}" CACHE STRING
    "ScoreP wrapper and flags")

  add_library(ScoreP::ScoreP INTERFACE IMPORTED)
  set_target_properties(ScoreP::ScoreP PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${ScoreP_INCLUDE_DIRS}")
endif()
