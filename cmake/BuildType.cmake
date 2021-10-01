# This file sets the default build type to Release, but allows the
# user to override using the usualy command line `-DCMAKE_BUILD_TYPE`
# option. Also respects multi-config generators like Visual Studio.
#
# Use like:
#
#     include(BuildType)
#
# Taken from https://blog.kitware.com/cmake-and-the-default-build-type/
# via the GS2 project https://bitbucket.org/gyrokinetics/utils/src/8.1-RC/cmake/BuildType.cmake

set(default_build_type "RelWithDebInfo")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()
