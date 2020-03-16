cmake_minimum_required(VERSION 3.13)

project(fci-wave-logn LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(fci-wave-logn fci-wave.cxx)
target_link_libraries(fci-wave-logn PRIVATE bout++::bout++)
