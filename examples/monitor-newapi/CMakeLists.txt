cmake_minimum_required(VERSION 3.13)

project(monitor-newapi LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(monitor-newapi monitor.cxx)
target_link_libraries(monitor-newapi PRIVATE bout++::bout++)
