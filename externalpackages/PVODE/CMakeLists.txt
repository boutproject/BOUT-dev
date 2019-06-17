cmake_minimum_required(VERSION 3.1...3.12)

if(${CMAKE_VERSION} VERSION_LESS 3.11)
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
  include/band.h
  include/cvdiag.h
  include/cvode.h
  include/cvspgmr.h
  include/iterativ.h
  include/llnlmath.h
  include/llnltyps.h
  include/nvector.h
  include/smalldense.h
  include/spgmr.h
  include/vector.h
  )

target_include_directories(pvode PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  )
target_link_libraries(pvode PUBLIC MPI::MPI_CXX)

add_library(pvpre
  precon/pvbbdpre.cpp
  precon/pvbbdpre.h
  precon/band.cpp
  precon/band.h
  )

target_include_directories(pvpre PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  )
target_link_libraries(pvpre PUBLIC MPI::MPI_CXX)

include(GNUInstallDirs)
install(TARGETS pvode pvpre
  EXPORT bout++Targets
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  )
