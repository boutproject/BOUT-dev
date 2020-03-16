cmake_minimum_required(VERSION 3.13)

project(blob2d-laplacexz LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(blob2d-laplacexz blob2d.cxx)
target_link_libraries(blob2d-laplacexz PRIVATE bout++::bout++)
