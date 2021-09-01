cmake_minimum_required(VERSION 3.13)

project(elmpb_outerloop LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

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


add_executable(elmpb_outerloop ${elmpb_sources})
target_link_libraries(elmpb_outerloop PRIVATE bout++::bout++)
if(BOUT_HAS_HYPRE)
   target_compile_definitions(elmpb_outerloop PUBLIC "BOUT_HAS_HYPRE")
   target_include_directories(elmpb_outerloop PUBLIC 
         ${HYPRE_DIR}/include 
         ${PROJECT_SOURCE_DIR}/include
         ${PROJECT_SOURCE_DIR}/externalpackages/fmt/include
         ${PROJECT_SOURCE_DIR}/externalpackages/mpark.variant/include
         )
   if(BOUT_USE_CUDA)
      list(APPEND elmpb_depends cuda)
      set_source_files_properties(${elmpb_sources} PROPERTIES LANGUAGE CUDA )
      message(STATUS "elmpb_outerloop needs CUDA LIBRARY DIR: ${BOUT_CUDA_LIB_DIR} and sets HYPRE_USING_CUDA") 
      target_compile_definitions(elmpb_outerloop PUBLIC "HYPRE_USING_CUDA;HYPRE_USING_UNIFIED_MEMORY") 
      target_include_directories(elmpb_outerloop PUBLIC ${BOUT_CUDA_LIB_DIR}/../include ) 
      target_compile_definitions(elmpb_outerloop PUBLIC "BOUT_USE_CUDA")
   endif ()
else()
   target_include_directories(elmpb_outerloop PUBLIC 
         ${PROJECT_SOURCE_DIR}/include
         ${PROJECT_SOURCE_DIR}/externalpackages/fmt/include
         ${PROJECT_SOURCE_DIR}/externalpackages/mpark.variant/include
         )
   if(BOUT_USE_CUDA)
      list(APPEND elmpb_depends cuda)
      set_source_files_properties(${elmpb_sources} PROPERTIES LANGUAGE CUDA )
      message(STATUS "elmpb_outerloop needs CUDA LIBRARY DIR: ${BOUT_CUDA_LIB_DIR}") 
      target_include_directories(elmpb_outerloop PUBLIC ${BOUT_CUDA_LIB_DIR}/../include ) 
      target_compile_definitions(elmpb_outerloop PUBLIC "BOUT_USE_CUDA")
   endif ()
endif ()

if(BOUT_HAS_RAJA)
   target_compile_definitions(elmpb_outerloop PUBLIC "BOUT_HAS_RAJA")
endif ()

message(STATUS "elmpb_depends: ${elmpb_depends}")

# Copy the documentation
bout_copy_file(README.md)

# Input files
bout_copy_file(data/BOUT.inp)
bout_copy_file(data-hypre/BOUT.inp)
bout_copy_file(data-nonlinear/BOUT.inp)
bout_copy_file(cbm18_dens8.grid_nx68ny64.nc)
bout_copy_file(cbm18_8_y064_x516_090309.nc)
