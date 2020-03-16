cmake_minimum_required(VERSION 3.13)

project(gyro-gem LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(gyro-gem gem.cxx)
target_link_libraries(gyro-gem PRIVATE bout++::bout++)
