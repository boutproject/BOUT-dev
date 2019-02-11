#ifndef __PAR_BNDRY_OP_H__
#define __PAR_BNDRY_OP_H__

#include "boundary_op.hxx"
#include "bout_types.hxx"
#include "parallel_boundary_region.hxx"
#include "utils.hxx"
#include "unused.hxx"

#include <utility>

//////////////////////////////////////////////////
// Base class

class BoundaryOpPar : public BoundaryOpBase {
public:
  BoundaryOpPar() : bndry(nullptr), real_value(0.), value_type(REAL) {}
  BoundaryOpPar(BoundaryRegionPar *region, std::shared_ptr<FieldGenerator> value)
      : bndry(region), gen_values(std::move(value)), value_type(GEN) {}
  BoundaryOpPar(BoundaryRegionPar *region, Field3D* value) :
    bndry(region),
    field_values(value),
    value_type(FIELD) {}
  BoundaryOpPar(BoundaryRegionPar *region, BoutReal value) :
    bndry(region),
    real_value(value),
    value_type(REAL) {}
  ~BoundaryOpPar() override {}

  // Note: All methods must implement clone, except for modifiers (see below)
  virtual BoundaryOpPar* clone(BoundaryRegionPar *UNUSED(region), const std::list<std::string> &UNUSED(args)) {return nullptr; }
  virtual BoundaryOpPar* clone(BoundaryRegionPar *UNUSED(region), Field3D *UNUSED(f)) {return nullptr; }

  virtual BoundaryOpPar*
  clone(BoundaryRegionPar* region, const std::list<std::string>& args,
        const std::map<std::string, std::string>& UNUSED(keywords)) {
    // If not implemented, call two-argument version
    return clone(region, args);
  }

  using BoundaryOpBase::apply;
  void apply(Field2D &UNUSED(f)) override {
    throw BoutException("Can't apply parallel boundary conditions to Field2D!");
  }
  void apply(Field2D &UNUSED(f), BoutReal UNUSED(t)) override {
    throw BoutException("Can't apply parallel boundary conditions to Field2D!");
  }

  BoundaryRegionPar *bndry;

protected:

  /// Possible ways to get boundary values
  std::shared_ptr<FieldGenerator>  gen_values;
  Field3D* field_values;
  BoutReal real_value;

  /// Where to take boundary values from - the generator, field or BoutReal
  enum ValueType { GEN, FIELD, REAL };
  const ValueType value_type;

  BoutReal getValue(int x, int y, int z, BoutReal t);
  BoutReal getValue(const BoundaryRegionPar &bndry, BoutReal t);

};

//////////////////////////////////////////////////
// Implementations

class BoundaryOpPar_dirichlet : public BoundaryOpPar {
public:
  BoundaryOpPar_dirichlet() : BoundaryOpPar(nullptr, 0.) {}
  BoundaryOpPar_dirichlet(BoundaryRegionPar *region) :
    BoundaryOpPar(region, 0.) {}
  BoundaryOpPar_dirichlet(BoundaryRegionPar *region, std::shared_ptr<FieldGenerator>  value) :
    BoundaryOpPar(region, value) {}
  BoundaryOpPar_dirichlet(BoundaryRegionPar *region, Field3D* value) :
    BoundaryOpPar(region, value) {}
  BoundaryOpPar_dirichlet(BoundaryRegionPar *region, BoutReal value) :
    BoundaryOpPar(region, value) {}
  BoundaryOpPar* clone(BoundaryRegionPar *region, const std::list<std::string> &args) override;
  BoundaryOpPar* clone(BoundaryRegionPar *region, Field3D *f) override;

  using BoundaryOpPar::apply;
  void apply(Field3D &f) override {return apply(f, 0);}
  void apply(Field3D &f, BoutReal t) override;

};

class BoundaryOpPar_dirichlet_O3 : public BoundaryOpPar {
public:
  BoundaryOpPar_dirichlet_O3() : BoundaryOpPar(nullptr, 0.) {}
  BoundaryOpPar_dirichlet_O3(BoundaryRegionPar *region) :
    BoundaryOpPar(region, 0.) {}
  BoundaryOpPar_dirichlet_O3(BoundaryRegionPar *region, std::shared_ptr<FieldGenerator>  value) :
    BoundaryOpPar(region, value) {}
  BoundaryOpPar_dirichlet_O3(BoundaryRegionPar *region, Field3D* value) :
    BoundaryOpPar(region, value) {}
  BoundaryOpPar_dirichlet_O3(BoundaryRegionPar *region, BoutReal value) :
    BoundaryOpPar(region, value) {}
  BoundaryOpPar* clone(BoundaryRegionPar *region, const std::list<std::string> &args) override;
  BoundaryOpPar* clone(BoundaryRegionPar *region, Field3D *f) override;

  using BoundaryOpPar::apply;
  void apply(Field3D &f) override {return apply(f, 0);}
  void apply(Field3D &f, BoutReal t) override;

};

class BoundaryOpPar_dirichlet_interp : public BoundaryOpPar {
public:
  BoundaryOpPar_dirichlet_interp() : BoundaryOpPar(nullptr, 0.) {}
  BoundaryOpPar_dirichlet_interp(BoundaryRegionPar *region) :
    BoundaryOpPar(region, 0.) {}
  BoundaryOpPar_dirichlet_interp(BoundaryRegionPar *region, std::shared_ptr<FieldGenerator>  value) :
    BoundaryOpPar(region, value) {}
  BoundaryOpPar_dirichlet_interp(BoundaryRegionPar *region, Field3D* value) :
    BoundaryOpPar(region, value) {}
  BoundaryOpPar_dirichlet_interp(BoundaryRegionPar *region, BoutReal value) :
    BoundaryOpPar(region, value) {}
  BoundaryOpPar* clone(BoundaryRegionPar *region, const std::list<std::string> &args) override;
  BoundaryOpPar* clone(BoundaryRegionPar *region, Field3D *f) override;

  using BoundaryOpPar::apply;
  void apply(Field3D &f) override {return apply(f, 0);}
  void apply(Field3D &f, BoutReal t) override;

};

class BoundaryOpPar_neumann : public BoundaryOpPar {
public:
  BoundaryOpPar_neumann() : BoundaryOpPar(nullptr, 0.) {}
  BoundaryOpPar_neumann(BoundaryRegionPar *region) :
    BoundaryOpPar(region, 0.) {}
  BoundaryOpPar_neumann(BoundaryRegionPar *region, std::shared_ptr<FieldGenerator>  value) :
    BoundaryOpPar(region, value) {}
  BoundaryOpPar_neumann(BoundaryRegionPar *region, Field3D* value) :
    BoundaryOpPar(region, value) {}
  BoundaryOpPar_neumann(BoundaryRegionPar *region, BoutReal value) :
    BoundaryOpPar(region, value) {}
  BoundaryOpPar* clone(BoundaryRegionPar *region, const std::list<std::string> &args) override;
  BoundaryOpPar* clone(BoundaryRegionPar *region, Field3D *f) override;

  using BoundaryOpPar::apply;
  void apply(Field3D &f) override {return apply(f, 0);}
  void apply(Field3D &f, BoutReal t) override;

};

#endif // __PAR_BNDRY_OP_H__
