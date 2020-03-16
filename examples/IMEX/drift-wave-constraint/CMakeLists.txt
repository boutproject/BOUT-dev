cmake_minimum_required(VERSION 3.13)

project(drift-wave-constraint LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(drift-wave-constraint test-drift.cxx)
target_link_libraries(drift-wave-constraint PRIVATE bout++::bout++)
