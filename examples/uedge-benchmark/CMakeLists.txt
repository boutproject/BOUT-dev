cmake_minimum_required(VERSION 3.13)

project(uedge-benchmark LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(uedge-benchmark ue_bmark.cxx)
target_link_libraries(uedge-benchmark PRIVATE bout++::bout++)
