cmake_minimum_required(VERSION 3.13)

project(hasegawa-wakatani LANGUAGES CXX)

if(NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(hasegawa-wakatani SOURCES hw.cxx)

if(BOUT_HAS_CUDA)
  set_source_files_properties(hw.cxx PROPERTIES LANGUAGE CUDA)
  target_link_options(hasegawa-wakatani
                      PRIVATE "-L${CUDAToolkit_ROOT}/targets/x86_64-linux/lib")
endif()
