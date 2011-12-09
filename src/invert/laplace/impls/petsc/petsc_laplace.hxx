
class LaplacePetsc;

#ifndef __PETSC_LAPLACE_H__
#define __PETSC_LAPLACE_H__

#include <invert_laplace.hxx>
#include <options.hxx>

class LaplacePetsc : public Laplacian {
public:
  LaplacePetsc(Options *opt = NULL) : Laplacian(opt) {}
  ~LaplacePetsc() {}
  
  void setCoefA(const Field2D &val) { A = val; }
  void setCoefC(const Field2D &val) { C = val; }
  void setCoefD(const Field2D &val) { D = val; }

  void setCoefA(const Field3D &val) { A = val; }
  void setCoefC(const Field3D &val) { C = val; }
  void setCoefD(const Field3D &val) { D = val; }
  
  const FieldPerp solve(const FieldPerp &b);
  const FieldPerp solve(const FieldPerp &b, const FieldPerp &x0);
  
private:
  Field3D A, C, D;
  
};

#endif __PETSC_LAPLACE_H__

