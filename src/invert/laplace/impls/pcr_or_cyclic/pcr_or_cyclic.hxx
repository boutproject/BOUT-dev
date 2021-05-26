#ifndef BOUT_PCR_CYCLIC_HXX
#define BOUT_PCR_CYCLIC_HXX

#include <invert_laplace.hxx>

#include <memory>

/// This is a wrapper around both the `LaplacePCR` _and_
/// `LaplaceCyclic` solvers: it defaults to using `LaplacePCR`, unless
/// one of the preconditions for using PCR aren't met, in which case
/// this solver falls back to using `LaplaceCyclic` instead
class LaplacePCRorCyclic : public Laplacian {
public:
  LaplacePCRorCyclic(Options* opt = nullptr, const CELL_LOC loc = CELL_CENTRE,
                     Mesh* mesh_in = nullptr);
  ~LaplacePCRorCyclic() = default;

  void setCoefA(const Field2D& val) override { pcr_or_cyclic->setCoefA(val); }
  void setCoefA(const Field3D& val) override { pcr_or_cyclic->setCoefA(val); }
  void setCoefA(BoutReal val) override { pcr_or_cyclic->setCoefA(val); }

  void setCoefC(const Field2D& val) override { pcr_or_cyclic->setCoefC(val); }
  void setCoefC(const Field3D& val) override { pcr_or_cyclic->setCoefC(val); }
  void setCoefC(BoutReal val) override { pcr_or_cyclic->setCoefC(val); }

  void setCoefC1(const Field2D& val) override { pcr_or_cyclic->setCoefC1(val); }
  void setCoefC1(const Field3D& val) override { pcr_or_cyclic->setCoefC1(val); }
  void setCoefC1(BoutReal val) override { pcr_or_cyclic->setCoefC1(val); }

  void setCoefC2(const Field2D& val) override { pcr_or_cyclic->setCoefC2(val); }
  void setCoefC2(const Field3D& val) override { pcr_or_cyclic->setCoefC2(val); }
  void setCoefC2(BoutReal val) override { pcr_or_cyclic->setCoefC2(val); }

  void setCoefD(const Field2D& val) override { pcr_or_cyclic->setCoefD(val); }
  void setCoefD(const Field3D& val) override { pcr_or_cyclic->setCoefD(val); }
  void setCoefD(BoutReal val) override { pcr_or_cyclic->setCoefD(val); }

  void setCoefEx(const Field2D& val) override { pcr_or_cyclic->setCoefEx(val); }
  void setCoefEx(const Field3D& val) override { pcr_or_cyclic->setCoefEx(val); }
  void setCoefEx(BoutReal val) override { pcr_or_cyclic->setCoefEx(val); }

  void setCoefEz(const Field2D& val) override { pcr_or_cyclic->setCoefEz(val); }
  void setCoefEz(const Field3D& val) override { pcr_or_cyclic->setCoefEz(val); }
  void setCoefEz(BoutReal val) override { pcr_or_cyclic->setCoefEz(val); }
  
  void setGlobalFlags(int f) override { pcr_or_cyclic->setGlobalFlags(f); }
  void setInnerBoundaryFlags(int f) override { pcr_or_cyclic->setInnerBoundaryFlags(f); }
  void setOuterBoundaryFlags(int f) override { pcr_or_cyclic->setOuterBoundaryFlags(f); }

  bool uses3DCoefs() const override { return pcr_or_cyclic->uses3DCoefs(); }

  FieldPerp solve(const FieldPerp& b) override { return pcr_or_cyclic->solve(b); }
  FieldPerp solve(const FieldPerp& rhs, const FieldPerp& x0) override {
    return pcr_or_cyclic->solve(rhs, x0);
  }

  Field3D solve(const Field3D& b) override { return pcr_or_cyclic->solve(b); }
  Field3D solve(const Field3D& rhs, const Field3D& x0) override {
    return pcr_or_cyclic->solve(rhs, x0);
  }

  Field2D solve(const Field2D& b) override { return pcr_or_cyclic->solve(b); }
  Field2D solve(const Field2D& rhs, const Field2D& x0) override {
    return pcr_or_cyclic->solve(rhs, x0);
  }

private:
  std::unique_ptr<Laplacian> pcr_or_cyclic;
};

namespace {
RegisterLaplace<LaplacePCRorCyclic> registerlaplacepcrorcyclic(LAPLACE_PCR_OR_CYCLIC);
}

#endif // BOUT_PCR_CYCLIC_HXX
