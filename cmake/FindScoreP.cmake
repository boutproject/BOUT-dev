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
#   ScoreP_LIBRARIES - List of libraries need to link against ScoreP
#   ScoreP_CXX_FLAGS - Compile definitions
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
#
# ----------------
# Part of this module (the bit that parses scorep-config for the libraries) was lifted from
# https://raw.githubusercontent.com/score-p/scorep_plugin_common/master/FindScorep.cmake
#
# Copyright (c) 2016, Technische Universität Dresden, Germany
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted
# provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions
#    and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions
#    and the following disclaimer in the documentation and/or other materials provided with the
#    distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse
#    or promote products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
# THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


find_program(ScoreP_CONFIG scorep-config)
mark_as_advanced(ScoreP_CONFIG)

get_filename_component(ScoreP_TMP "${ScoreP_CONFIG}" DIRECTORY)
get_filename_component(ScoreP_EXEC_LOCATION "${ScoreP_TMP}" DIRECTORY)

if (ScoreP_DEBUG)
  message(STATUS "[ ${CMAKE_CURRENT_LIST_FILE}:${CMAKE_CURRENT_LIST_LINE} ] "
    " ScoreP_CONFIG = ${ScoreP_CONFIG}"
    " ScoreP_EXEC_LOCATION = ${ScoreP_EXEC_LOCATION}")
endif()

if(ScoreP_CONFIG)
  message(STATUS "SCOREP library found. (using ${ScoreP_CONFIG})")

  execute_process(COMMAND ${ScoreP_CONFIG} "--user" "--nocompiler" "--cppflags"
    OUTPUT_VARIABLE ScoreP_CONFIG_FLAGS)

  string(REGEX MATCHALL "-I[^ ]*" ScoreP_CONFIG_INCLUDES "${ScoreP_CONFIG_FLAGS}")
  foreach(inc ${ScoreP_CONFIG_INCLUDES})
    string(SUBSTRING ${inc} 2 -1 inc)
    list(APPEND ScoreP_INCLUDE_DIRS ${inc})
  endforeach()

  if (ScoreP_DEBUG)
    message(STATUS "[ ${CMAKE_CURRENT_LIST_FILE}:${CMAKE_CURRENT_LIST_LINE} ] "
      " ScoreP_INCLUDE_DIRS = ${ScoreP_INCLUDE_DIRS}")
  endif()

  string(REGEX MATCHALL "(^| +)-[^I][^ ]*" ScoreP_CONFIG_CXXFLAGS "${ScoreP_CONFIG_FLAGS}")
  foreach(flag ${ScoreP_CONFIG_CXXFLAGS})
    string(STRIP ${flag} flag)
    list(APPEND ScoreP_CXX_FLAGS ${flag})
  endforeach()

  if (ScoreP_DEBUG)
    message(STATUS "[ ${CMAKE_CURRENT_LIST_FILE}:${CMAKE_CURRENT_LIST_LINE} ] "
      " ScoreP_CXX_FLAGS = ${ScoreP_CXX_FLAGS}")
  endif()

  unset(ScoreP_CONFIG_FLAGS)
  unset(ScoreP_CONFIG_INCLUDES)
  unset(ScoreP_CONFIG_CXXFLAGS)

  execute_process(COMMAND ${ScoreP_CONFIG} "--user" "--nocompiler" "--ldflags"
    OUTPUT_VARIABLE _LINK_LD_ARGS)
  string( REPLACE " " ";" _LINK_LD_ARGS ${_LINK_LD_ARGS} )
  foreach( _ARG ${_LINK_LD_ARGS} )
    if(${_ARG} MATCHES "^-L")
      STRING(REGEX REPLACE "^-L" "" _ARG ${_ARG})
      SET(ScoreP_LINK_DIRS ${ScoreP_LINK_DIRS} ${_ARG})
    endif()
  endforeach()

  if (ScoreP_DEBUG)
    message(STATUS "[ ${CMAKE_CURRENT_LIST_FILE}:${CMAKE_CURRENT_LIST_LINE} ] "
      " ScoreP_LINK_DIRS = ${ScoreP_LINK_DIRS}")
  endif()

  execute_process(COMMAND ${ScoreP_CONFIG} "--user" "--nocompiler" "--libs"
    OUTPUT_VARIABLE _LINK_LD_ARGS)
  string( REPLACE " " ";" _LINK_LD_ARGS ${_LINK_LD_ARGS} )
  foreach( _ARG ${_LINK_LD_ARGS} )
    if(${_ARG} MATCHES "^-l")
      string(REGEX REPLACE "^-l" "" _ARG ${_ARG})
      find_library(_SCOREP_LIB_FROM_ARG NAMES ${_ARG}
        PATHS
        ${ScoreP_LINK_DIRS}
        )
      if(_SCOREP_LIB_FROM_ARG)
        set(ScoreP_LIBRARIES ${ScoreP_LIBRARIES} ${_SCOREP_LIB_FROM_ARG})
      endif()
      unset(_SCOREP_LIB_FROM_ARG CACHE)
    endif()
  endforeach()

  if (ScoreP_DEBUG)
    message(STATUS "[ ${CMAKE_CURRENT_LIST_FILE}:${CMAKE_CURRENT_LIST_LINE} ] "
      " ScoreP_LIBRARIES = ${ScoreP_LIBRARIES}")
  endif()

endif()

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args(ScoreP DEFAULT_MSG
    ScoreP_CONFIG
    ScoreP_LIBRARIES
    ScoreP_INCLUDE_DIRS
  )

if (ScoreP_FOUND AND NOT TARGET ScoreP::ScoreP)
  add_library(ScoreP::ScoreP UNKNOWN IMPORTED)
  set_target_properties(ScoreP::ScoreP PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${ScoreP_INCLUDE_DIRS}"
    IMPORTED_LINK_INTERFACE_LIBRARIES "${ScoreP_LIBRARIES}"
    INTERFACE_INCLUDE_DEFINITIONS "${ScoreP_CXX_FLAGS}"
    IMPORTED_LINK_INTERFACE_LANGUAGES "C"
    )
endif()
