#include "gtest/gtest.h"

#include "bout/conservative_flux_div.hxx"

// Simple test to verify that conservative_flux_div compiles and links correctly
TEST(ConservativeFluxDivCompileTest, HeaderIncludesCompile) {
  // This test just verifies that we can include the header
  // and that the implementation exists in the library
  SUCCEED();
}

// Test that the namespace and function signatures are correct
TEST(ConservativeFluxDivCompileTest, FunctionSignatures) {
  // Verify the function exists in the FV namespace
  // We can't call it without a proper mesh, but we can verify it exists
  using FuncType1 = Field3D (*)(const FaceField3D&, bool);
  using FuncType2 = Field3D (*)(const FaceField3D&, const Field3D*, const Field3D*, const Field3D*);
  
  // These lines verify the function signatures are correct
  FuncType1 func1 = &FV::ConservativeFluxDiv;
  FuncType2 func2 = &FV::ConservativeFluxDiv;
  
  // Avoid unused variable warnings
  (void)func1;
  (void)func2;
  
  SUCCEED();
}