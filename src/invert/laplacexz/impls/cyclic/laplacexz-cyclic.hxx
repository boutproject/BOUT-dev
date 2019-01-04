#ifndef BOUT_LAPLACEXZ_CYCLIC_H
#define BOUT_LAPLACEXZ_CYCLIC_H

#include "bout/invert/laplacexz.hxx"
#include "cyclic_reduction.hxx"
#include "dcomplex.hxx"
#include <map>
#include <memory>
#include <vector>

#include "bout/array.hxx"
#include "bout_types.hxx"
#include "field2d.hxx"
#include "field3d.hxx"
#include "utils.hxx"

class Mesh;
class Options;

class LaplaceXZcyclic : public LaplaceXZ {
public:
  LaplaceXZcyclic(Mesh *m = nullptr, Options *options = nullptr,
      const CELL_LOC loc = CELL_CENTRE);
  ~LaplaceXZcyclic() {}
  
  using LaplaceXZ::setCoefs;
  void setCoefs(const Field2D &A, const Field2D &B) override;
  
  using LaplaceXZ::solve;
  Field3D solve(const Field3D &b, const Field3D &x0) override;
private:
  int xstart, xend;
  int nmode, nloc, nsys;
  Matrix<dcomplex> acoef, bcoef, ccoef, xcmplx, rhscmplx;
  Array<dcomplex> k1d, k1d_2;
  std::unique_ptr<CyclicReduce<dcomplex>> cr; ///< Tridiagonal solver

  int inner_boundary_flags; ///< Flags to set inner boundary condition
  int outer_boundary_flags; ///< Flags to set outer boundary condition
};

#endif // BOUT_LAPLACEXZ_CYCLIC_H
