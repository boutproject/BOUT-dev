cmake_minimum_required(VERSION 3.13)

project(monitor LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(monitor monitor.cxx)
target_link_libraries(monitor PRIVATE bout++::bout++)
