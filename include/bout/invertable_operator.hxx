/**************************************************************************
 * Invert arbitrary linear global operation using PETSc.
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

namespace bout {
namespace inversion {
template <typename T>
class InvertableOperator;
};
};

#ifndef __INVERTABLE_OPERATOR_H__
#define __INVERTABLE_OPERATOR_H__

#ifdef BOUT_HAS_PETSC

#include "bout/traits.hxx"
#include <bout/mesh.hxx>
#include <bout/sys/timer.hxx>
#include <boutcomm.hxx>
#include <boutexception.hxx>
#include <globals.hxx>
#include <msg_stack.hxx>
#include <options.hxx>
#include <output.hxx>

#include <petscksp.h>

#include <bout/petsclib.hxx>

#endif

namespace bout {
namespace inversion {

#ifdef BOUT_HAS_PETSC

/// No-op function to use as a default -- may wish to remove once testing phase complete
template <typename T>
T identity(const T& in) {
  AUTO_TRACE();
  return in;
};

/// Pack a PetscVec from a Field<T>
template <typename T>
PetscErrorCode fieldToPetscVec(const T& in, Vec out) {
  TRACE("fieldToPetscVec<T>");
  Timer timer("invertable_operator_packing");

  PetscScalar* vecData;

  auto ierr = VecGetArray(out, &vecData);
  CHKERRQ(ierr);

  int counter = 0;

  // Should explore ability to OpenMP this
  BOUT_FOR_SERIAL(i, in.getRegion("RGN_WITHBNDRIES")) {
    vecData[counter] = in[i];
    counter++;
  }

  ierr = VecRestoreArray(out, &vecData);
  CHKERRQ(ierr);

  return ierr;
}

/// Pack a Field<T> from a PetscVec
template <typename T>
PetscErrorCode petscVecToField(Vec in, T& out) {
  TRACE("petscVecToField<T>");
  Timer timer("invertable_operator_packing");

  const PetscScalar* vecData;

  auto ierr = VecGetArrayRead(in, &vecData);
  CHKERRQ(ierr);

  int counter = 0;

  // Should explore ability to OpenMP this
  BOUT_FOR_SERIAL(i, out.getRegion("RGN_WITHBNDRIES")) {
    out[i] = vecData[counter];
    counter++;
  }

  ierr = VecRestoreArrayRead(in, &vecData);
  CHKERRQ(ierr);

  return ierr;
};

/// Class to define an invertable operator. Provides interface to PETSc routines
/// for solving A.x = b
template <typename T>
class InvertableOperator {
  static_assert(
      bout::utils::is_Field<T>::value,
      "InvertableOperator must be templated with one of FieldPerp, Field2D or Field3D");

public:
  /// What type of field does the operator take?
  using data_type = T;

  /// The signature of the functor that applies the operator.
  using function_signature = std::function<T(const T&)>;

  /// Almost empty constructor -- currently don't actually use Options for anything
  InvertableOperator(const function_signature& func = identity<T>,
                     Options* optIn = nullptr, Mesh* localmeshIn = nullptr)
      : operatorFunction(func), preconditionerFunction(func),
        opt(optIn == nullptr ? optIn
                             : Options::getRoot()->getSection("invertableOperator")),
        localmesh(localmeshIn == nullptr ? bout::globals::mesh : localmeshIn),
        lib(opt) {
    AUTO_TRACE();
  };

  /// Destructor just has to cleanup the PETSc owned objects.
  ~InvertableOperator() {
    TRACE("InvertableOperator<T>::destructor");

    KSPDestroy(&ksp);
    MatDestroy(&matOperator);
    MatDestroy(&matPreconditioner);
    VecDestroy(&rhs);
    VecDestroy(&lhs);
  };

  /// Allow the user to override the existing function
  /// Note by default we set the preconditioner function to match this
  /// as this is the usual mode of operation. If the user doesn't want to
  /// do this they can set alsoSetPreconditioner to false.
  void setOperatorFunction(const function_signature& func,
                           bool alsoSetPreconditioner = true) {
    TRACE("InvertableOperator<T>::setOperatorFunction");
    operatorFunction = func;
    if (alsoSetPreconditioner) {
      preconditionerFunction = func;
    }
  }

  /// Allow the user to override the existing preconditioner function
  void setPreconditionerFunction(const function_signature& func) {
    TRACE("InvertableOperator<T>::setPreconditionerFunction");
    preconditionerFunction = func;
  }

  /// Provide a way to apply the operator to a Field
  T operator()(const T& input) {
    TRACE("InvertableOperator<T>::operator()");
    return operatorFunction(input);
  }

  /// Provide a synonym for applying the operator to a Field
  T apply(const T& input) {
    AUTO_TRACE();
    return operator()(input);
  }

  /// Sets up the PETSc objects required for inverting the operator
  /// Currently also takes the functor that applies the operator this class
  /// represents. Not actually required by any of the setup so this should
  /// probably be moved to a separate place (maybe the constructor).
  PetscErrorCode setup() {
    TRACE("InvertableOperator<T>::setup");

    Timer timer("invertable_operator_setup");
    if (doneSetup) {
      throw BoutException(
          "Trying to call setup on an InvertableOperator instance that has "
          "already been setup.");
    }

    // Add the RGN_WITHBNDRIES region to the mesh. Requires RGN_NOBNDRY to be defined.
    if (std::is_same<Field3D, T>::value) {
      if (not localmesh->hasRegion3D("RGN_WITHBNDRIES")) {
        // This avoids all guard cells and corners but includes boundaries
        // Note we probably don't want to include periodic boundaries as these
        // are essentially just duplicate points so should be careful here (particularly
        // in y)
        // to only include unique points
        Region<Ind3D> nocorner3D = localmesh->getRegion3D("RGN_NOBNDRY");
        if (!localmesh->periodicX) {
          if (localmesh->firstX())
            nocorner3D += Region<Ind3D>(0, localmesh->xstart - 1, localmesh->ystart,
                                        localmesh->yend, 0, localmesh->LocalNz - 1,
                                        localmesh->LocalNy, localmesh->LocalNz,
                                        localmesh->maxregionblocksize);
          if (localmesh->lastX())
            nocorner3D += Region<Ind3D>(
                localmesh->LocalNx - localmesh->xstart, localmesh->LocalNx - 1,
                localmesh->ystart, localmesh->yend, 0, localmesh->LocalNz - 1,
                localmesh->LocalNy, localmesh->LocalNz, localmesh->maxregionblocksize);
        }
        if (localmesh->firstY() or localmesh->lastY()) {
          for (int ix = localmesh->xstart; ix <= localmesh->xend; ix++) {
            if (not localmesh->periodicY(ix)) {
              if (localmesh->firstY())
                nocorner3D +=
                    Region<Ind3D>(ix, ix, 0, localmesh->ystart - 1, 0,
                                  localmesh->LocalNz - 1, localmesh->LocalNy,
                                  localmesh->LocalNz, localmesh->maxregionblocksize);
              if (localmesh->lastY())
                nocorner3D += Region<Ind3D>(
                    ix, ix, localmesh->LocalNy - localmesh->ystart,
                    localmesh->LocalNy - 1, 0, localmesh->LocalNz - 1, localmesh->LocalNy,
                    localmesh->LocalNz, localmesh->maxregionblocksize);
            }
          }
        }

        nocorner3D.unique();
        localmesh->addRegion3D("RGN_WITHBNDRIES", nocorner3D);
      }

    } else if (std::is_same<Field2D, T>::value) {
      if (not localmesh->hasRegion2D("RGN_WITHBNDRIES")) {
        // This avoids all guard cells and corners but includes boundaries
        Region<Ind2D> nocorner2D = localmesh->getRegion2D("RGN_NOBNDRY");
        if (!localmesh->periodicX) {
          if (localmesh->firstX())
            nocorner2D += Region<Ind2D>(0, localmesh->xstart - 1, localmesh->ystart,
                                        localmesh->yend, 0, 0, localmesh->LocalNy, 1,
                                        localmesh->maxregionblocksize);
          if (localmesh->lastX())
            nocorner2D +=
                Region<Ind2D>(localmesh->LocalNx - localmesh->xstart,
                              localmesh->LocalNx - 1, localmesh->ystart, localmesh->yend,
                              0, 0, localmesh->LocalNy, 1, localmesh->maxregionblocksize);
        }
        if (localmesh->firstY() or localmesh->lastY()) {
          for (int ix = localmesh->xstart; ix <= localmesh->xend; ix++) {
            if (not localmesh->periodicY(ix)) {
              if (localmesh->firstY())
                nocorner2D +=
                    Region<Ind2D>(ix, ix, 0, localmesh->ystart - 1, 0, 0,
                                  localmesh->LocalNy, 1, localmesh->maxregionblocksize);
              if (localmesh->lastY())
                nocorner2D +=
                    Region<Ind2D>(ix, ix, localmesh->LocalNy - localmesh->ystart,
                                  localmesh->LocalNy - 1, 0, 0, localmesh->LocalNy, 1,
                                  localmesh->maxregionblocksize);
            }
          }
        }
        nocorner2D.unique();
        localmesh->addRegion2D("RGN_WITHBNDRIES", nocorner2D);
      }

    } else if (std::is_same<FieldPerp, T>::value) {
      if (not localmesh->hasRegionPerp("RGN_WITHBNDRIES")) {
        // This avoids all guard cells and corners but includes boundaries
        Region<IndPerp> nocornerPerp = localmesh->getRegionPerp("RGN_NOBNDRY");
        if (!localmesh->periodicX) {
          if (localmesh->firstX())
            nocornerPerp +=
                Region<IndPerp>(0, localmesh->xstart - 1, 0, 0, 0, localmesh->LocalNz - 1,
                                1, localmesh->LocalNz, localmesh->maxregionblocksize);
          if (localmesh->lastX())
            nocornerPerp +=
                Region<IndPerp>(localmesh->LocalNx - localmesh->xstart,
                                localmesh->LocalNx - 1, 0, 0, 0, localmesh->LocalNz - 1,
                                1, localmesh->LocalNz, localmesh->maxregionblocksize);
        }
        nocornerPerp.unique();
        localmesh->addRegionPerp("RGN_WITHBNDRIES", nocornerPerp);
      }

    } else {
      throw BoutException("Invalid template type provided to InvertableOperator");
    }

    // Hacky way to determine the local size for now
    PetscInt nlocal = 0;
    {
      T tmp(localmesh);
      nlocal = tmp.getRegion("RGN_WITHBNDRIES").size();
    }

    PetscInt nglobal = PETSC_DETERMINE; // Depends on type of T

    /// Create the shell matrix representing the operator to invert
    /// Note we currently pass "this" as the Matrix context
    PetscInt ierr = MatCreateShell(BoutComm::get(), nlocal, nlocal, nglobal, nglobal,
                                   this, &matOperator);
    CHKERRQ(ierr);

/// Create vectors compatible with matrix
#if PETSC_VERSION_LT(3, 6, 0)
    ierr = MatGetVecs(matOperator, &rhs, &lhs);
#else
    ierr = MatCreateVecs(matOperator, &rhs, &lhs);
#endif
    CHKERRQ(ierr);

    /// Zero out the lhs vector as values used for initial guess
    /// in solve step. Don't need to initialise rhs as this is
    /// done in the callback function using the passed Field.
    ierr = VecSet(lhs, 0.0);
    CHKERRQ(ierr);

    /// Now register Matrix_multiply operation
    ierr = MatShellSetOperation(matOperator, MATOP_MULT, (void (*)())(functionWrapper));
    CHKERRQ(ierr);

    /// Create the shell matrix representing the operator to invert
    /// Note we currently pass "this" as the Matrix context
    ierr = MatCreateShell(BoutComm::get(), nlocal, nlocal, nglobal, nglobal, this,
                          &matPreconditioner);
    CHKERRQ(ierr);

    /// Now register Matrix_multiply operation
    ierr = MatShellSetOperation(matPreconditioner, MATOP_MULT,
                                (void (*)())(preconditionerWrapper));
    CHKERRQ(ierr);

    /// Now create and setup the linear solver with the matrix
    lib.createKSPWithOptions(BoutComm::get(), &ksp);

#if PETSC_VERSION_LT(3, 5, 0)
    /// Need to provide a MatStructure flag in versions <3.5. This details if we expect
    /// the preconditioner matrix structure to vary between calls to KSPSolve.
    /// Safest but slowest option is DIFFERENT_NONZERO_PATTERN but can probably usually
    /// use SAME_PRECONDITIONER.
    ierr =
        KSPSetOperators(ksp, matOperator, matPreconditioner, DIFFERENT_NONZERO_PATTERN);
#else
    ierr = KSPSetOperators(ksp, matOperator, matPreconditioner);
#endif
    CHKERRQ(ierr);

    /// By default allow a non-zero initial guess as this is probably the
    /// most helpful mode of operation. To disable this user can pass
    /// `-invertable_ksp_initial_guess_nonzero false` on the command line or simply
    /// use `result = operator.invert(rhs, T{})` which will lead to an initial
    /// guess of whatever is in the instance of T that is passed.
    ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    CHKERRQ(ierr);

    /// Allow options to be set on command line using a --invertable_ksp_* prefix.
    ierr = KSPSetOptionsPrefix(ksp, "invertable_");
    CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp);
    CHKERRQ(ierr);

    /// Do required setup so solve can proceed in invert
    ierr = KSPSetUp(ksp);
    CHKERRQ(ierr);

    doneSetup = true;

    return ierr;
  };

  // Pack values into lhs vec before solving.
  // This may be the way to set the initial guess
  // but suspect it's not as there are KSPGuess objects
  // to deal with.
  T invert(const T& rhsField, const T& guess) {
    AUTO_TRACE();
    auto ierr = fieldToPetscVec(guess, lhs);
    CHKERRQ(ierr);
    return invert(rhsField);
  }

  /// Triggers the solve of A.x = b for x, where b = rhs and A is the matrix
  /// representation
  /// of the operator we represent. Should probably provide an overload or similar as a
  /// way of setting the initial guess.
  T invert(const T& rhsField) {
    TRACE("InvertableOperator<T>::invert");
    Timer timer("invertable_operator_invert");

    if (!doneSetup) {
      throw BoutException("Trying to call invert on an InvertableOperator instance that "
                          "has not been setup.");
    }

    ASSERT2(localmesh == rhsField.getMesh());

    // rhsField to rhs
    auto ierr = fieldToPetscVec(rhsField, rhs);
    CHKERRQ(ierr);

    /// Do the solve with solution stored in lhs
    /// Note: the values in lhs on input are used as the initial guess
    /// provided KSPSetInitialGuessNonzero has been called as true (if
    /// not then lhs is zeroed on entry). As the results in lhs persist
    /// between calls to invert (as a class member rather than local scope)
    /// we automatically provide the previous solution as the initial guess
    /// for subsequent solves.
    ierr = KSPSolve(ksp, rhs, lhs);
    CHKERRQ(ierr);

    KSPConvergedReason reason;
    ierr = KSPGetConvergedReason(ksp, &reason);
    if (reason <= 0) {
      throw BoutException("KSPSolve failed with reason %d.", reason);
    }

    // Probably want to remove the following in the long run
    output_debug << "KSPSolve finished with converged reason : " << reason << endl;

    // lhs to lhsField -- first make the output field and ensure it has space allocated
    T lhsField{emptyFrom(rhsField)};

    ierr = petscVecToField(lhs, lhsField);
    CHKERRQ(ierr);

    return lhsField;
  };

  /// With checks enabled provides a convience routine to check that
  /// applying the registered function on the calculated inverse gives
  /// back the initial values.
  bool verify(const T& rhsIn, BoutReal tol = 1.0e-5) {
    TRACE("InvertableOperator<T>::verify");

    T result = invert(rhsIn);
    localmesh->communicate(result);
    const T applied = operator()(result);
    const BoutReal maxDiff = max(abs(applied - rhsIn), true);
    if (maxDiff >= tol) {
      output_debug << "Maximum difference in verify is " << maxDiff << endl;
      output_debug << "Max rhs is " << max(abs(rhsIn), true) << endl;
      output_debug << "Max applied is " << max(abs(applied), true) << endl;
      output_debug << "Max result is " << max(abs(result), true) << endl;
    };
    return maxDiff < tol;
  };

  /// Reports the time spent in various parts of InvertableOperator. Note
  /// that as the Timer "labels" are not unique to an instance the time
  /// reported is summed across all different instances.
  static void reportTime() {
    TRACE("InvertableOperator<T>::reportTime");
    BoutReal time_setup = Timer::resetTime("invertable_operator_setup");
    BoutReal time_invert = Timer::resetTime("invertable_operator_invert");
    BoutReal time_packing = Timer::resetTime("invertable_operator_packing");

    BoutReal time_operate = Timer::resetTime("invertable_operator_operate");
    output_warn << "InvertableOperator timing :: Setup " << time_setup;
    output_warn << " , Invert(packing, operation) " << time_invert << "(";
    output_warn << time_packing << ", ";
    output_warn << time_operate << "). Total : " << time_setup + time_invert << endl;
  };

private:
  // PETSc objects
  Mat matOperator, matPreconditioner;
  Vec rhs, lhs;
  KSP ksp;

  /// The function that represents the operator that we wish to invert
  function_signature operatorFunction = identity<T>;

  /// The function that represents the preconditioner for the operator that we wish to
  /// invert
  function_signature preconditionerFunction = identity<T>;

  // Internal types
  Options* opt = nullptr;    // Do we need this?
  Mesh* localmesh = nullptr; //< To ensure we can create T on the right mesh
  bool doneSetup = false;

  // To ensure PETSc has been setup -- a bit noisy if creating/destroying
  // InvertableOperator,
  // maybe this should be static to avoid this but then how do we initialise it?
  PetscLib lib; // Do we need this?

  /// Wrapper that gets a pointer to the parent InvertableOperator instance
  /// from the Matrix m and uses this to get the actual function to call.
  /// Copies data from v1 into a field of type T, calls the function on this and then
  /// copies the result into the v2 argument.
  static PetscErrorCode functionWrapper(Mat m, Vec v1, Vec v2) {
    TRACE("InvertableOperator<T>::functionWrapper");
    InvertableOperator<T>* ctx;
    auto ierr = MatShellGetContext(m, &ctx);
    T tmpField(ctx->localmesh);
    tmpField.allocate();
    ierr = petscVecToField(v1, tmpField);
    CHKERRQ(ierr);
    // Need following communicate if operator() uses guard cells, i.e. differential
    // operator. Could delegate to the user function but then need to remove const
    // from signature of the function (function_signature) likely involving a copy.
    // @TODO : Consider removing the communicate and introduce requirement for user
    // function to communicate if required. This would be neater as currently result
    // of invert needs explicitly communicating if we want to apply the operator to
    // it, for example (e.g. see verify).
    ctx->localmesh->communicate(tmpField);
    T tmpField2 = ctx->operator()(tmpField);
    // This communicate is required in case operator() ends up not setting
    // all periodic boundaries correctly (possibly -- need to check?)
    // @TODO : Consider need for this communicate. Could communicate at the
    // end of the user routine.
    ctx->localmesh->communicate(tmpField2);
    ierr = fieldToPetscVec(tmpField2, v2);
    CHKERRQ(ierr);
    return ierr;
  }

  /// Wrapper that gets a pointer to the parent InvertableOperator instance
  /// from the Matrix m and uses this to get the actual function to call.
  /// Copies data from v1 into a field of type T, calls the function on this and then
  /// copies the result into the v2 argument.
  static PetscErrorCode preconditionerWrapper(Mat m, Vec v1, Vec v2) {
    TRACE("InvertableOperator<T>::functionWrapper");
    InvertableOperator<T>* ctx;
    auto ierr = MatShellGetContext(m, &ctx);
    T tmpField(ctx->localmesh);
    tmpField.allocate();
    ierr = petscVecToField(v1, tmpField);
    CHKERRQ(ierr);
    // Need following communicate if operator() uses guard cells, i.e. differential
    // operator. Could delegate to the user function but then need to remove const
    // from signature of the function (function_signature) likely involving a copy.
    // @TODO : Consider removing the communicate and introduce requirement for user
    // function to communicate if required. This would be neater as currently result
    // of invert needs explicitly communicating if we want to apply the operator to
    // it, for example (e.g. see verify).
    ctx->localmesh->communicate(tmpField);
    T tmpField2 = ctx->preconditionerFunction(tmpField);
    // This communicate is required in case operator() ends up not setting
    // all periodic boundaries correctly (possibly -- need to check?)
    // @TODO : Consider need for this communicate. Could communicate at the
    // end of the user routine.
    ctx->localmesh->communicate(tmpField2);
    ierr = fieldToPetscVec(tmpField2, v2);
    CHKERRQ(ierr);
    return ierr;
  }
};

#else

template <typename T>
class InvertableOperator {
public:
};

#endif // PETSC
};
};

#endif // HEADER GUARD
