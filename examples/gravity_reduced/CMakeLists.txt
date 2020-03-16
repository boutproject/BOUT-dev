cmake_minimum_required(VERSION 3.13)

project(gravity_reduced LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(gravity_reduced gravity_reduced.cxx)
target_link_libraries(gravity_reduced PRIVATE bout++::bout++)
