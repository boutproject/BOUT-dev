
#include <bout/invert/laplacexz.hxx>
#include <cyclic_reduction.hxx>
#include <dcomplex.hxx>

class LaplaceXZcyclic : public LaplaceXZ {
public:
  LaplaceXZcyclic(Mesh *m, Options *options);
  ~LaplaceXZcyclic();
  
  void setCoefs(const Field2D &A, const Field2D &B);
  
  Field3D solve(const Field3D &b, const Field3D &x0);
private:
  Mesh *mesh;   ///< The mesh this operates on, provides metrics and communication

  int xstart, xend;
  int nmode, nloc, nsys;
  dcomplex **acoef, **bcoef, **ccoef, **xcmplx, **rhscmplx, *k1d;
  CyclicReduce<dcomplex> *cr; ///< Tridiagonal solver
};
