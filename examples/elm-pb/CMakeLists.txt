cmake_minimum_required(VERSION 3.13)


project(elmpb LANGUAGES CXX C)

set(elmpb_sources
   elm_pb.cxx
)

if(BOUT_HAS_HYPRE)
   set(elmpb_depends bout++ PETSc HYPRE)
else()
   set(elmpb_depends bout++ PETSc)
endif ()

if(BOUT_HAS_RAJA)
   list(APPEND elmpb_depends RAJA)
   find_package(RAJA REQUIRED)
endif ()

if(BOUT_HAS_UMPIRE)
   list(APPEND elmpb_depends umpire)
   find_package(umpire REQUIRED)
endif ()

add_executable(elmpb ${elmpb_sources})
target_link_libraries(elmpb PRIVATE bout++::bout++)
if(BOUT_HAS_HYPRE)
   target_compile_definitions(elmpb PUBLIC "BOUT_HAS_HYPRE")
   target_include_directories(elmpb PUBLIC 
         ${HYPRE_DIR}/include 
         ${PROJECT_SOURCE_DIR}/include
         ${PROJECT_SOURCE_DIR}/externalpackages/fmt/include
         ${PROJECT_SOURCE_DIR}/externalpackages/mpark.variant/include
         )
   if(BOUT_USE_CUDA)
      list(APPEND elmpb_depends cuda)
      set_source_files_properties(${elmpb_sources} PROPERTIES LANGUAGE CUDA )
      message(STATUS "elmpb needs CUDA LIBRARY DIR: ${BOUT_CUDA_LIB_DIR} and sets HYPRE_USING_CUDA") 
      target_compile_definitions(elmpb PUBLIC "HYPRE_USING_CUDA;HYPRE_USING_UNIFIED_MEMORY") 
      target_include_directories(elmpb PUBLIC ${BOUT_CUDA_LIB_DIR}/../include ) 
      target_compile_definitions(elmpb PUBLIC "BOUT_USE_CUDA")
   endif ()
else()
   target_include_directories(elmpb PUBLIC 
         ${PROJECT_SOURCE_DIR}/include
         ${PROJECT_SOURCE_DIR}/externalpackages/fmt/include
         ${PROJECT_SOURCE_DIR}/externalpackages/mpark.variant/include
         )
   if(BOUT_USE_CUDA)
      list(APPEND elmpb_depends cuda)
      set_source_files_properties(${elmpb_sources} PROPERTIES LANGUAGE CUDA )
      message(STATUS "elmpb needs CUDA LIBRARY DIR: ${BOUT_CUDA_LIB_DIR}") 
      target_include_directories(elmpb PUBLIC ${BOUT_CUDA_LIB_DIR}/../include ) 
      target_compile_definitions(elmpb PUBLIC "BOUT_USE_CUDA")
   endif ()
endif ()

if(BOUT_HAS_RAJA)
   target_compile_definitions(elmpb PUBLIC "BOUT_HAS_RAJA")
endif ()
message(STATUS "elmpb_depends: ${elmpb_depends}")
