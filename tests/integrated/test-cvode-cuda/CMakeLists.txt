bout_add_integrated_test(
  test-cvode-cuda
  SOURCES test_cvode_cuda.cxx
  REQUIRES BOUT_HAS_CVODE BOUT_HAS_SUNDIALS_CUDA
  USE_DATA_BOUT_INP USE_RUNTEST
)
if(TARGET test-cvode-cuda)
  set_source_files_properties(test_cvode_cuda.cxx PROPERTIES LANGUAGE CUDA)
  set_target_properties(test-cvode-cuda PROPERTIES CUDA_STANDARD 20)
  set_tests_properties(
    test-cvode-cuda PROPERTIES LABELS "GPU;CUDA" PROCESSORS 2
  )
endif()
