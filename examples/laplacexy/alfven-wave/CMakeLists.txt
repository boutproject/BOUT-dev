cmake_minimum_required(VERSION 3.13)

project(laplacexy-alfven-wave LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(laplacexy-alfven-wave alfven.cxx)
target_link_libraries(laplacexy-alfven-wave PRIVATE bout++::bout++)
