cmake_minimum_required(VERSION 3.13)

project(elm_pb LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(elm_pb elm_pb.cxx)
target_link_libraries(elm_pb PRIVATE bout++::bout++)
