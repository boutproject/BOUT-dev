cmake_minimum_required(VERSION 3.13)

project(subsampling LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(subsampling monitor.cxx)
target_link_libraries(subsampling PRIVATE bout++::bout++)
