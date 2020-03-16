cmake_minimum_required(VERSION 3.13)

project(constraints-alfven-wave LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(constraints-alfven-wave alfven.cxx)
target_link_libraries(constraints-alfven-wave PRIVATE bout++::bout++)
