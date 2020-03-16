cmake_minimum_required(VERSION 3.13)

project(lapd-drift LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(lapd-drift lapd_drift.cxx)
target_link_libraries(lapd-drift PRIVATE bout++::bout++)
