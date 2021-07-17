cmake_minimum_required(VERSION 3.13)

if (${CMAKE_MINOR_VERSION} VERSION_GREATER_EQUAL "18") 
  cmake_policy(SET CMP0104 OLD) # we set cuda arch explicitly otherwise need CUDA_ARCHITECTURES for each target
endif ()


project(hw3d LANGUAGES CXX C)

set(hw_sources
   hw.cxx
)

set(hw_depends bout++ PETSc)
if(BOUT_HAS_RAJA)
   list(APPEND hw_depends RAJA)
   find_package(RAJA REQUIRED)
endif ()

if(BOUT_HAS_UMPIRE)
   list(APPEND hw_depends umpire)
   find_package(umpire REQUIRED)
endif ()

add_executable(hw ${hw_sources})
target_link_libraries(hw PRIVATE bout++::bout++)
target_include_directories(hw PUBLIC 
      ${PROJECT_SOURCE_DIR}/include
      ${PROJECT_SOURCE_DIR}/externalpackages/fmt/include
      ${PROJECT_SOURCE_DIR}/externalpackages/mpark.variant/include
      )
if(BOUT_USE_CUDA)
   list(APPEND hw_depends cuda)
   set_source_files_properties(${hw_sources} PROPERTIES LANGUAGE CUDA )
   message(STATUS "hw needs CUDA LIBRARY DIR: ${BOUT_CUDA_LIB_DIR}") 
   target_include_directories(hw PUBLIC ${BOUT_CUDA_LIB_DIR}/../include ) 
   target_compile_definitions(hw PUBLIC "BOUT_USE_CUDA")
endif ()

if(BOUT_HAS_RAJA)
   target_compile_definitions(hw PUBLIC "BOUT_HAS_RAJA")
endif ()

message(STATUS "hw_depends: ${hw_depends}")
