cmake_minimum_required(VERSION 3.13)

project(test_laplacexy LANGUAGES CXX)

find_package(bout++ REQUIRED)
add_executable(test_laplacexy test-laplacexy.cxx)
target_link_libraries(test_laplacexy PRIVATE bout++::bout++)
