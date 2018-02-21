class Laplace3D;

#ifndef __LAPLACE3D_H__
#define __LAPLACE3D_H__

///////////////////////////////////////////////////////////////////
// Interface to 3D Laplacian solver
// 
// Solves the following:
// 
// A*x + B * Delp^2(x) + Div( C * Delp(x)) + D * Del^2(x) = rhs
//
// where Del^2 is the full 3D Laplacian
//       Delp^2 is the full perpendicular Laplacian
//              including Y derivatives
//
///////////////////////////////////////////////////////////////////

#include <bout/options.hxx>
#include <bout/field2d.hxx>
#include <bout/field3d.hxx>

class Laplace3D {
public:
  Laplace3D(Options *opt = NULL) : flags(0) {}
  virtual ~Laplace3D() {}
  
  virtual void setCoefA(const Field2D &f) = 0;
  virtual void setCoefA(const Field3D &f) {setCoefA(DC(f));}
  virtual void setCoefA(BoutReal f) {setCoefA(Field2D(f));}
  
  virtual void setCoefB(const Field2D &f) = 0;
  virtual void setCoefB(const Field3D &f) {setCoefB(DC(f));}
  virtual void setCoefB(BoutReal f) {setCoefB(Field2D(f));}
  
  virtual void setCoefC(const Field2D &f) = 0;
  virtual void setCoefC(const Field3D &f) {setCoefC(DC(f));}
  virtual void setCoefC(BoutReal f) {setCoefC(Field2D(f));}
  
  virtual void setCoefD(const Field2D &f) = 0;
  virtual void setCoefD(const Field3D &f) {setCoefD(DC(f));}
  virtual void setCoefD(BoutReal f) {setCoefD(Field2D(f));}

  virtual void setFlags(int f) {flags = f;}
  
  virtual const Field3D solve(const Field3D &b) = 0;
  virtual const Field2D solve(const Field2D &b) {return solve(DC(Field3D(b)));}
  
  virtual const Field3D solve(const Field3D &b, const Field3D &x0) { return solve(b); }
  virtual const Field2D solve(const Field2D &b, const Field2D &x0) { return solve(b); }
  
  static Laplace3D* create(Options *opt = NULL);
protected:
  int flags;
  
};

#endif //__LAPLACE3D_H__

