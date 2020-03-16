cmake_minimum_required(VERSION 3.13)

project(finite-volume-test LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(finite-volume-test finite_volume.cxx)
target_link_libraries(finite-volume-test PRIVATE bout++::bout++)
