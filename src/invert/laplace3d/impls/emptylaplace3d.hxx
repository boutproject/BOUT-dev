
#ifndef __EMPTYLAPLACE3D_H__
#define __EMPTYLAPLACE3D_H__

#include <bout/invert/laplace3d.hxx>
#include <boutexception.hxx>

class EmptyLaplace3D : public Laplace3D {
public:
  EmptyLaplace3D(Options *opt = nullptr) {
    throw BoutException("Laplace3D solver not available");
  }
  void setCoefA(const Field2D &f) {}
  void setCoefB(const Field2D &f) {}
  void setCoefC(const Field2D &f) {}
  void setCoefD(const Field2D &f) {}
 
  const Field3D solve(const Field3D &b) {return b;}
};

#endif // __EMPTYLAPLACE3D_H__
