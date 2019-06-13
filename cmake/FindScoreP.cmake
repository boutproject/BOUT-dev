find_program(SCOREP_EXECUTABLE scorep)
mark_as_advanced(SCOREP_EXECUTABLE)

if (NOT SCOREP_EXECUTABLE)
  message(FATAL_ERROR "Score-P requested but not executable not found.
Please supply the path using -DCMAKE_PROGRAM_PATH=/path/to/scorep/bin/directory")
endif()

get_filename_component(SCOREP_TMP "${SCOREP_EXECUTABLE}" DIRECTORY)
get_filename_component(SCOREP_EXEC_LOCATION "${SCOREP_TMP}" DIRECTORY)

find_path(ScoreP_INCLUDE_DIR
  NAMES scorep/SCOREP_User.h
  DOC "Score-P include directory"
  HINTS "${SCOREP_EXEC_LOCATION}"
  PATH_SUFFIXES "include"
  )
mark_as_advanced(ScoreP_INCLUDE_DIR)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ScoreP
  REQUIRED_VARS ScoreP_INCLUDE_DIR
  )

if (ScoreP_FOUND AND NOT TARGET ScoreP::ScoreP)
  set(ScoreP_FLAGS "--user;--nocompiler" CACHE STRING "ScoreP wrapper flags")
  set(ScoreP_COMPILER_LAUNCHER "${SCOREP_EXECUTABLE};${ScoreP_FLAGS}" CACHE STRING
    "ScoreP wrapper and flags")

  add_library(ScoreP::ScoreP INTERFACE IMPORTED)
  set_target_properties(ScoreP::ScoreP PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${ScoreP_INCLUDE_DIR}")
endif()
