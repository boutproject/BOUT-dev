cmake_minimum_required(VERSION 3.13)

project(staggered_grid LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(staggered_grid test_staggered.cxx)
target_link_libraries(staggered_grid PRIVATE bout++::bout++)
