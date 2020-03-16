cmake_minimum_required(VERSION 3.13)

project(dalf3 LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(dalf3 dalf3.cxx)
target_link_libraries(dalf3 PRIVATE bout++::bout++)
