cmake_minimum_required(VERSION 3.9..3.12)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
else()
  cmake_policy(VERSION 3.12)
endif()

project(PVODE
  DESCRIPTION "ODE Solver"
  VERSION 0.1
  LANGUAGES CXX)

find_package(MPI REQUIRED)

add_library(pvode
  source/cvode.cpp
  source/nvector.cpp
  source/llnlmath.cpp
  source/cvspgmr.cpp
  source/spgmr.cpp
  source/iterativ.cpp
  source/cvdiag.cpp
  source/smalldense.cpp
  include/pvode/band.h
  include/pvode/cvdiag.h
  include/pvode/cvode.h
  include/pvode/cvspgmr.h
  include/pvode/iterativ.h
  include/pvode/llnlmath.h
  include/pvode/llnltyps.h
  include/pvode/nvector.h
  include/pvode/smalldense.h
  include/pvode/spgmr.h
  include/pvode/vector.h
  )

target_include_directories(pvode PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/pvode>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  )
target_link_libraries(pvode PUBLIC MPI::MPI_CXX)

add_library(pvpre
  include/pvode/pvbbdpre.h
  precon/pvbbdpre.cpp
  precon/band.cpp
  precon/band.h
  )

target_include_directories(pvpre PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/pvode>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  )
target_link_libraries(pvpre PUBLIC MPI::MPI_CXX)

include(GNUInstallDirs)
install(TARGETS pvode pvpre
  EXPORT PVODETargets
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  )
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  PVODEConfigVersion.cmake
  VERSION ${PACKAGE_VERSION}
  COMPATIBILITY SameMajorVersion
  )

install(EXPORT PVODETargets
  FILE PVODETargets.cmake
  NAMESPACE PVODE::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/PVODE"
  )

export(EXPORT PVODETargets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/PVODEConfig.cmake"
  NAMESPACE PVODE::
  )

export(PACKAGE PVODE)
