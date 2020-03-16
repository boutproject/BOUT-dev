cmake_minimum_required(VERSION 3.13)

project(finite-volume-fluid LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(finite-volume-fluid fluid.cxx)
target_link_libraries(finite-volume-fluid PRIVATE bout++::bout++)
