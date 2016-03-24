/************************************************************
 * Perpendicular Laplacian inversion using PETSc
 *
 *
 ************************************************************/

class LaplaceXZpetsc;

#ifndef __LAPLACEXZ_PETSC_H__
#define __LAPLACEXZ_PETSC_H__

#include <bout/invert/laplacexz.hxx>

#ifndef BOUT_HAS_PETSC

#include <boutexception.hxx>
class LaplaceXZpetsc : public LaplaceXZ {
public:
  LaplaceXZpetsc(Mesh *m, Options *options) : LaplaceXZ(m, options) {
    throw BoutException("No PETSc LaplaceXY solver available");
  }
  void setCoefs(const Field2D &A, const Field2D &B) {}
  Field3D solve(const Field3D &b, const Field3D &x0) {}
private:
};

#else // BOUT_HAS_PETSC

#include <bout/petsclib.hxx>

class LaplaceXZpetsc : public LaplaceXZ {
public:
  /*!
   * Constructor
   */
  LaplaceXZpetsc(Mesh *m, Options *options);

  /*!
   * Destructor
   */
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
  vector<YSlice> slice;

  Vec xs, bs;        ///< Solution and RHS vectors

  Mesh *mesh;   ///< The mesh this operates on, provides metrics and communication

  int reuse_limit; ///< How many times can the preconditioner be reused?
  int reuse_count; ///< How many times has it been reused?

  bool coefs_set; ///< Have coefficients been set?

  static const int INVERT_AC_GRAD  = 2;  // Use zero neumann (NOTE: AC is a misnomer)
  static const int INVERT_SET      = 16; // Set boundary to x0 value
  static const int INVERT_RHS      = 32; // Set boundary to b value
  #ifdef CHECK
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
