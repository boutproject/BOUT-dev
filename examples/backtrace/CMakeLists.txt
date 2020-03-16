cmake_minimum_required(VERSION 3.13)

project(backtrace LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(backtrace backtrace.cxx)
target_link_libraries(backtrace PRIVATE bout++::bout++)
