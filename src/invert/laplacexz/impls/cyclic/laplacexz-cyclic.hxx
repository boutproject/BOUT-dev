#ifndef LAPLACEXZ_CYCLIC_HXX
#define LAPLACEXZ_CYCLIC_HXX

#include "bout/build_defines.hxx"
#include "bout/invert/laplacexz.hxx"

#if BOUT_USE_METRIC_3D

namespace {
RegisterUnavailableLaplaceXZ registerlaplacexzcylic{
    "cyclic", "BOUT++ was configured with 3D metrics"};
} // namespace

#else

#include "bout/utils.hxx"
#include <bout/cyclic_reduction.hxx>
#include <bout/dcomplex.hxx>
#include <bout/globals.hxx>

class LaplaceXZcyclic : public LaplaceXZ {
public:
  LaplaceXZcyclic(Mesh* m = nullptr, Options* options = nullptr,
                  const CELL_LOC loc = CELL_CENTRE);
  ~LaplaceXZcyclic() {}

  using LaplaceXZ::setCoefs;
  void setCoefs(const Field2D& A, const Field2D& B) override;

  using LaplaceXZ::solve;
  Field3D solve(const Field3D& b, const Field3D& x0) override;

private:
  int xstart, xend;
  int nmode, nloc, nsys;
  Matrix<dcomplex> acoef, bcoef, ccoef, xcmplx, rhscmplx;
  Array<dcomplex> k1d, k1d_2;
  std::unique_ptr<CyclicReduce<dcomplex>> cr; ///< Tridiagonal solver

  int inner_boundary_flags; ///< Flags to set inner boundary condition
  int outer_boundary_flags; ///< Flags to set outer boundary condition
};

namespace {
RegisterLaplaceXZ<LaplaceXZcyclic> registerlaplacexzcylic{"cyclic"};
} // namespace

#endif // BOUT_USE_METRIC_3D

#endif // LAPLACEXZ_CYCLIC_HXX
