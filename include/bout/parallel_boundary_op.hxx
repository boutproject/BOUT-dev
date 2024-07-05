#ifndef BOUT_PAR_BNDRY_OP_H
#define BOUT_PAR_BNDRY_OP_H

#include "bout/boundary_op.hxx"
#include "bout/bout_types.hxx"
#include "bout/field_factory.hxx"
#include "bout/parallel_boundary_region.hxx"
#include "bout/unused.hxx"
#include "bout/utils.hxx"

#include <utility>

//////////////////////////////////////////////////
// Base class

class BoundaryOpPar : public BoundaryOpBase {
public:
  BoundaryOpPar() = default;
  BoundaryOpPar(BoundaryRegionPar* region, std::shared_ptr<FieldGenerator> value)
      : bndry(region), gen_values(std::move(value)), value_type(ValueType::GEN) {}
  BoundaryOpPar(BoundaryRegionPar* region, Field3D* value)
      : bndry(region), field_values(value), value_type(ValueType::FIELD) {}
  BoundaryOpPar(BoundaryRegionPar* region, BoutReal value)
      : bndry(region), real_value(value), value_type(ValueType::REAL) {}
  BoundaryOpPar(BoundaryRegionPar* region)
      : bndry(region), real_value(0.), value_type(ValueType::REAL) {}
  ~BoundaryOpPar() override = default;

  // Note: All methods must implement clone, except for modifiers (see below)
  virtual BoundaryOpPar* clone(BoundaryRegionPar* region,
                               const std::list<std::string>& args) = 0;
  virtual BoundaryOpPar* clone(BoundaryRegionPar* region, Field3D* f) = 0;
  virtual BoundaryOpPar*
  clone(BoundaryRegionPar* region, const std::list<std::string>& args,
        const std::map<std::string, std::string>& UNUSED(keywords)) {
    // If not implemented, call two-argument version
    return clone(region, args);
  }

  BoundaryRegionPar* bndry{nullptr};

protected:
  /// Possible ways to get boundary values
  std::shared_ptr<FieldGenerator> gen_values;
  Field3D* field_values{nullptr};
  BoutReal real_value{0.};

  /// Where to take boundary values from - the generator, field or BoutReal
  enum class ValueType { GEN, FIELD, REAL };
  const ValueType value_type{ValueType::REAL};

  BoutReal getValue(const BoundaryRegionPar& bndry, BoutReal t);
};

template <class T, bool isNeumann = false>
class BoundaryOpParTemp : public BoundaryOpPar {
public:
  using BoundaryOpPar::BoundaryOpPar;

  using BoundaryOpPar::clone;

  // Note: All methods must implement clone, except for modifiers (see below)
  BoundaryOpPar* clone(BoundaryRegionPar* region,
                       const std::list<std::string>& args) override {
    if (!args.empty()) {
      try {
        real_value = stringToReal(args.front());
        return new T(region, real_value);
      } catch (const BoutException&) {
        std::shared_ptr<FieldGenerator> newgen = nullptr;
        // First argument should be an expression
        newgen = FieldFactory::get()->parse(args.front());
        return new T(region, newgen);
      }
    }

    return new T(region);
  }

  BoundaryOpPar* clone(BoundaryRegionPar* region, Field3D* f) override {
    return new T(region, f);
  }

  using BoundaryOpBase::apply;
  void apply(Field2D& UNUSED(f)) final {
    throw BoutException("Can't apply parallel boundary conditions to Field2D!");
  }
  void apply(Field2D& UNUSED(f), BoutReal UNUSED(t)) final {
    throw BoutException("Can't apply parallel boundary conditions to Field2D!");
  }
  void apply(Field3D& f) override { return apply(f, 0); }

  void apply(Field3D& f, BoutReal t) override {
    f.ynext(bndry->dir).allocate(); // Ensure unique before modifying

    auto dy = f.getCoordinates()->dy;

    for (bndry->first(); !bndry->isDone(); bndry->next()) {
      BoutReal value = getValue(*bndry, t);
      if (isNeumann) {
        value *= dy[bndry->ind()];
      }
      static_cast<T*>(this)->apply_stencil(f, bndry, value);
    }
  }
};

//////////////////////////////////////////////////
// Implementations

class BoundaryOpPar_dirichlet_o1 : public BoundaryOpParTemp<BoundaryOpPar_dirichlet_o1> {
public:
  using BoundaryOpParTemp::BoundaryOpParTemp;
  static void apply_stencil(Field3D& f, const BoundaryRegionPar* bndry, BoutReal value) {
    bndry->dirichlet_o1(f, value);
  }
};

class BoundaryOpPar_dirichlet_o2 : public BoundaryOpParTemp<BoundaryOpPar_dirichlet_o2> {
public:
  using BoundaryOpParTemp::BoundaryOpParTemp;
  static void apply_stencil(Field3D& f, const BoundaryRegionPar* bndry, BoutReal value) {
    bndry->dirichlet_o2(f, value);
  }
};

class BoundaryOpPar_dirichlet_o3 : public BoundaryOpParTemp<BoundaryOpPar_dirichlet_o3> {
public:
  using BoundaryOpParTemp::BoundaryOpParTemp;
  static void apply_stencil(Field3D& f, const BoundaryRegionPar* bndry, BoutReal value) {
    bndry->dirichlet_o3(f, value);
  }
};

class BoundaryOpPar_neumann_o1
    : public BoundaryOpParTemp<BoundaryOpPar_neumann_o1, true> {
public:
  using BoundaryOpParTemp::BoundaryOpParTemp;
  static void apply_stencil(Field3D& f, const BoundaryRegionPar* bndry, BoutReal value) {
    bndry->neumann_o1(f, value);
  }
};

class BoundaryOpPar_neumann_o2
    : public BoundaryOpParTemp<BoundaryOpPar_neumann_o2, true> {
public:
  using BoundaryOpParTemp::BoundaryOpParTemp;
  static void apply_stencil(Field3D& f, const BoundaryRegionPar* bndry, BoutReal value) {
    bndry->neumann_o2(f, value);
  }
};

class BoundaryOpPar_neumann_o3
    : public BoundaryOpParTemp<BoundaryOpPar_neumann_o3, true> {
public:
  using BoundaryOpParTemp::BoundaryOpParTemp;
  static void apply_stencil(Field3D& f, const BoundaryRegionPar* bndry, BoutReal value) {
    bndry->neumann_o3(f, value);
  }
};

#endif // BOUT_PAR_BNDRY_OP_H
