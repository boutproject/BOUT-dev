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
  LaplaceXZpetsc(Mesh *m, Options *options) {
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
  
  void setCoefs(const Field2D &A, const Field2D &B);

  /*!
   * Solve Laplacian in X-Z
   */
  Field3D solve(const Field3D &b, const Field3D &x0);
  
private:
  PetscLib lib;     ///< Requires PETSc library
  Mat MatA;         ///< Matrix to be inverted
  Vec xs, bs;       ///< Solution and RHS vectors
  KSP ksp;          ///< Krylov Subspace solver
  PC pc;            ///< Preconditioner

  Mesh *mesh;   ///< The mesh this operates on, provides metrics and communication

  /*!
   * Number of grid points on this processor
   */
  int localSize();

  /*!
   * Return the communicator
   */
  MPI_Comm communicator();
};


#endif // BOUT_HAS_PETSC
#endif // __LAPLACEXZ_PETSC_H__
