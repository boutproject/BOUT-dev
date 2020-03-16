cmake_minimum_required(VERSION 3.13)

project(orszag-tang LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(orszag-tang mhd.cxx)
target_link_libraries(orszag-tang PRIVATE bout++::bout++)
