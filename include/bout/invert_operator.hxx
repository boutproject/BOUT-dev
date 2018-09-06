/**************************************************************************
 * Invert arbitrary global operatorn using PETSc. to invert
 *
 **************************************************************************
 * Copyright 2018 D. Dickinson
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
 *
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/

template<typename T> class InvertOperator;

#ifndef __INVERT_OPERATOR_H__
#define __INVERT_OPERATOR_H__

#ifdef BOUT_HAS_PETSC

#include <petscksp.h>

#include <bout/petsclib.hxx>

#include <boutcomm.hxx>
#include <boutexception.hxx>
#include <globals.hxx>
#include <options.hxx>
#include <output.hxx>

/// No-op function to use as a default -- may wish to remove once testing phase complete
template<typename T>
T identity(const T& in){return in;};

template<typename T>
class InvertOperator {
  static_assert(std::is_base_of<Field3D, T>::value || std::is_base_of<Field2D, T>::value || std::is_base_of<FieldPerp, T>::value,
                "InvertOperator must be templated with one of FieldPerp, Field2D or Field3D"); 
 
 public:

  using data_type = T; // What type of field
  using function_signature = std::function<T(const T&)>; // The required signature of the callback function representing the operator
  using function_cast = void(*)(void); // Signature that we have to cast the function to for PETSc
    
  InvertOperator(Options *opt = nullptr, Mesh *localmesh  = nullptr) : opt(opt ? opt : Options::getRoot()->getSection("invertOperator")),
						  localmesh(localmesh ? localmesh : mesh), doneSetup(false) {};

  ~InvertOperator(){
    
#if CHECK > 1
    output_info<<endl;
    output_info<<"Destroying KSP object in InvertOperator with properties: "<<endl;
    KSPView(ksp, PETSC_VIEWER_STDOUT_SELF);
    output_info<<endl;    
#endif
    
    KSPDestroy(&ksp);
    MatDestroy(&matOperator);
    VecDestroy(&rhs);
    VecDestroy(&lhs);
  };
  
  PetscErrorCode setup(function_signature funcIn = identity<T>);
  T invert(const T &rhs);

  /// With checks enabled provides a convience routine to check that
  /// applying the registered function on the calculated inverse gives
  /// back the initial values.
  bool verify(const T &rhs, BoutReal tol = 1.0e-5){
    TRACE("InvertOperator<T>::verify");
#if CHECK > 1
    const T result = invert(rhs);
    const T applied = func(result);
    const BoutReal maxDiff = max(abs(applied-rhs),true);
#if CHECK > 2
    if (maxDiff >= tol) {
      output_info << "Maximum difference in verify is "<<maxDiff<<endl;
      output_info << "Max rhs is "<<max(abs(rhs),true)<<endl;
      output_info << "Max applied is "<<max(abs(applied),true)<<endl;
      output_info << "Max result is "<<max(abs(result),true)<<endl;
    };
#endif
    return maxDiff < tol;
#else
    return true;
#endif    
  };

  /// Wrapper that gets a pointer to the parent InvertOperator instance
  /// from the Matrix m and uses this to get the actual function to call.
  /// Copies data from v1 into a field of type T, calls the function on this and then
  /// copies the result into the v2 argument.
  static PetscErrorCode functionWrapper(Mat m, Vec v1, Vec v2){
    InvertOperator<T> *ctx;
    auto ierr = MatShellGetContext(m, &ctx);
    T tmpField(ctx->localmesh); tmpField.allocate();
    petscVecToField(v1, tmpField);
    T tmpField2 = ctx->func(tmpField);
    fieldToPetscVec(tmpField2, v2);
    return ierr;
  }
  
  /// The function that represents the operator that we wish to invert
  function_signature func;
  
  
 private:
  // PETSc types
  Mat matOperator;
  Vec rhs, lhs;
  KSP ksp;

  // Internal types
  Options *opt; //Do we need this?
  Mesh *localmesh; //< To ensure we can create T on the right mesh
  bool doneSetup = false;
  
  // To ensure PETSc has been setup
  PetscLib lib; //Do we need this?
};

template<typename T>
PetscErrorCode InvertOperator<T>::setup(function_signature funcIn) {
  TRACE("InvertOperator<T>::setup");  
  if (doneSetup) {
    throw BoutException("Trying to call setup on an InvertOperator instance that has already been setup.");
  }

  func = funcIn;
  
  PetscInt ierr;
  
  // Hacky way to determine the local size for now
  PetscInt nlocal = 0;
  {
    T tmp(localmesh);
    /// @TODO : Replace this with a BOUT_FOR type loop or alternative (T.size()?)
    for(const auto &i: tmp){ nlocal++;}
  }
  
  PetscInt nglobal = PETSC_DETERMINE; // Depends on type of T

  /// Create the shell matrix representing the operator to invert
  /// Note we currently pass "this" as the Matrix context
  ierr = MatCreateShell(BoutComm::get(), nlocal, nlocal, nglobal, nglobal, this, &matOperator);
  CHKERRQ(ierr); // Can't call this in constructor as includes a return statement

  /// Create vectors compatible with matrix
  ierr = MatCreateVecs(matOperator, &rhs, &lhs); // Older versions may need to use MatGetVecs
  CHKERRQ(ierr);
  
  /// Now register Matrix_multiply operation
  //  ierr = MatShellSetOperation(matOperator, MATOP_MULT, (function_cast)(identityFunc));
  ierr = MatShellSetOperation(matOperator, MATOP_MULT, (function_cast)(functionWrapper));
  CHKERRQ(ierr);

  /// Now create and setup the linear solver with the matrix
  ierr = KSPCreate(BoutComm::get(), &ksp);
  CHKERRQ(ierr);    
  ierr = KSPSetOperators(ksp, matOperator, matOperator);
  CHKERRQ(ierr);

  /// Allow options to be set on command line using a --invert_ksp_* prefix.
  ierr = KSPSetOptionsPrefix(ksp,"invert_");
  CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);
  CHKERRQ(ierr);

  /// Do required setup so solve can proceed in invert
  ierr = KSPSetUp(ksp);
  CHKERRQ(ierr);

  doneSetup = true;
};

template<typename T>
T InvertOperator<T>::invert(const T& rhsField){
  TRACE("InvertOperator<T>::invert");
  
  if (!doneSetup) {
    throw BoutException("Trying to call invert on an InvertOperator instance that has not been setup.");
  }

  ASSERT2(localmesh == rhsField.getMesh());

  // rhsField to rhs
  fieldToPetscVec(rhsField, rhs);
  
  /// Do the solve with solution stored in lhs 
  auto ierr = KSPSolve(ksp, rhs, lhs);
  CHKERRQ(ierr);

  KSPConvergedReason reason;
  ierr = KSPGetConvergedReason( ksp, &reason );
  if(reason <=0 ){
    throw BoutException("KSPSolve failed with reason %d.",reason);
  }
  
#if CHECK > 3
  output_info << "KSPSolve finished with converged reason : "<<reason<<endl;
#endif
  
  // lhs to lhsField -- first make the output field and ensure it has space allocated
  T lhsField(localmesh);
  lhsField.allocate();

  petscVecToField(lhs, lhsField);

  return lhsField;
}

template<typename T>
PetscErrorCode fieldToPetscVec(const T& in, Vec out){
  TRACE("fieldToPetscVec<T>");
  
  PetscScalar *vecData;

  auto ierr = VecGetArray(out, &vecData);
  CHKERRQ(ierr);

  int counter = 0;

  /// @TODO : Replace this with a BOUT_FOR type loop
  for(const auto &i: in){
    vecData[counter] = in[i];
    counter++;
  }
  
  ierr = VecRestoreArray(out, &vecData);
  CHKERRQ(ierr);  
}

template<typename T>
PetscErrorCode petscVecToField(Vec in, T& out){
  TRACE("petscVecToField<T>");
  
  const PetscScalar *vecData;

  auto ierr = VecGetArrayRead(in, &vecData);
  CHKERRQ(ierr);

  int counter = 0;

  /// @TODO : Replace this with a BOUT_FOR type loop  
  for(const auto &i: out){
    out[i] = vecData[counter];
    counter++;
  }
  
  ierr = VecRestoreArrayRead(in, &vecData);
  CHKERRQ(ierr);  
}

#endif // PETSC
#endif // HEADER GUARD
