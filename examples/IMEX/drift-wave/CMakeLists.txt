cmake_minimum_required(VERSION 3.13)

project(drift-wave LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(test-drift test-drift.cxx)
target_link_libraries(test-drift PRIVATE bout++::bout++)
