cmake_minimum_required(VERSION 3.13)

project(laplacexy-simple LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(laplacexy-simple test-laplacexy.cxx)
target_link_libraries(laplacexy-simple PRIVATE bout++::bout++)
