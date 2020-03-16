cmake_minimum_required(VERSION 3.13)

project(tokamak-2fluid LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(tokamak-2fluid 2fluid.cxx)
target_link_libraries(tokamak-2fluid PRIVATE bout++::bout++)
