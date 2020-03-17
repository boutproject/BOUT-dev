cmake_minimum_required(VERSION 3.13)

project(blob2d LANGUAGES CXX)

find_package(bout++ REQUIRED)
add_executable(blob2d blob2d.cxx)
target_link_libraries(blob2d PRIVATE bout++::bout++)
