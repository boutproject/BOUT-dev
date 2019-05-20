/*
 * 3D Laplacian solver using PETSc
 */

#ifndef BOUT_HAS_PETSC
// No PETSc available, so define as an empty laplace3d

#include "../emptylaplace3d.hxx"
using Laplace3DPetsc = EmptyLaplace3D;

#else // BOUT_HAS_PETSC

class Laplace3DPetsc;

#ifndef __PETSC_LAPLACE3D_H__
#define __PETSC_LAPLACE3D_H__

#include <bout/invert/laplace3d.hxx>

class Laplace3DPetsc : public Laplace3D {
public:
  Laplace3DPetsc(Options *opt = nullptr);
  ~Laplace3DPetsc();
  
  void setCoefA(const Field2D &f) {A = f;}
  void setCoefB(const Field2D &f) {B = f;}
  void setCoefC(const Field2D &f) {C = f;}
  void setCoefD(const Field2D &f) {D = f;}
  
  void setCoefA(const Field3D &f) {A = f;}
  void setCoefB(const Field3D &f) {B = f;}
  void setCoefC(const Field3D &f) {C = f;}
  void setCoefD(const Field3D &f) {D = f;}
  
  const Field3D solve(const Field3D &b);
private:
  Field3D A,B,C,D; // Coefficients
};

#endif // __PETSC_LAPLACE3D_H__

#endif // BOUT_HAS_PETSC

