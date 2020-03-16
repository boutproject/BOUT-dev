cmake_minimum_required(VERSION 3.13)

project(eigen-box LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(eigen-box eigen-box.cxx)
target_link_libraries(eigen-box PRIVATE bout++::bout++)
