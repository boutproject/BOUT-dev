/************************************************************
 * Perpendicular Laplacian inversion using PETSc
 *
 *
 ************************************************************/

class LaplaceXZpetsc;

#ifndef __LAPLACEXZ_PETSC_H__
#define __LAPLACEXZ_PETSC_H__

#include "bout/build_config.hxx"
#include "bout/invert/laplacexz.hxx"

#if not BOUT_HAS_PETSC

namespace {
RegisterUnavailableLaplaceXZ
    registerlaplacexzpetsc("petsc", "BOUT++ was not configured with PETSc");
}

#else

#include <bout/petsclib.hxx>

namespace {
RegisterLaplaceXZ<LaplaceXZpetsc> registerlaplacexzpetsc{"petsc"};
}

class LaplaceXZpetsc : public LaplaceXZ {
public:
  LaplaceXZpetsc(Mesh *m = nullptr, Options *options = nullptr, const CELL_LOC loc = CELL_CENTRE);
  ~LaplaceXZpetsc();

  void setCoefs(const Field3D &A, const Field3D &B);

  void setCoefs(const Field2D &A, const Field2D &B) {
    setCoefs(Field3D(A), Field3D(B));
  }

  /*!
   * Solve Laplacian in X-Z
   */
  Field3D solve(const Field3D &b, const Field3D &x0);

private:
  PetscLib lib;      ///< Requires PETSc library

  /*!
   * Data for a single Y slice
   */
  struct YSlice {
    int yindex; /// Y index
    Mat MatA;  ///< Matrix to be inverted
    Mat MatP;  ///< Matrix for preconditioner
    KSP ksp;   ///< Krylov Subspace solver context
  };
  std::vector<YSlice> slice;

  Vec xs, bs;        ///< Solution and RHS vectors

  int reuse_limit; ///< How many times can the preconditioner be reused?
  int reuse_count; ///< How many times has it been reused?

  bool coefs_set; ///< Have coefficients been set?
  
  #if CHECK > 0
    // Currently implemented flags
    static const int implemented_boundary_flags =   INVERT_AC_GRAD
                                                  + INVERT_SET
                                                  + INVERT_RHS
                                                  ;
  #endif

  int inner_boundary_flags; ///< Flags to set inner boundary condition
  int outer_boundary_flags; ///< Flags to set outer boundary condition
};


#endif // BOUT_HAS_PETSC
#endif // __LAPLACEXZ_PETSC_H__
