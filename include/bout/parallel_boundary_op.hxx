#ifndef BOUT_PAR_BNDRY_OP_H
#define BOUT_PAR_BNDRY_OP_H

#include "bout/boundary_op.hxx"
#include "bout/boundary_region_iter.hxx"
#include "bout/bout_types.hxx"
#include "bout/boutexception.hxx"
#include "bout/field3d.hxx"
#include "bout/field_factory.hxx"
#include "bout/sys/expressionparser.hxx"
#include "bout/unused.hxx"

#include <memory>

//////////////////////////////////////////////////
// Base class

class BoundaryOpPar : public BoundaryOpBase {
public:
  BoundaryOpPar() = default;
  BoundaryOpPar(bout::boundary::BoundaryRegionFCI* region,
                std::shared_ptr<FieldGenerator> value)
      : bndry(region), gen_values(std::move(value)), value_type(ValueType::GEN) {}
  BoundaryOpPar(bout::boundary::BoundaryRegionFCI* region, Field3D* value)
      : bndry(region), field_values(value), value_type(ValueType::FIELD) {}
  BoundaryOpPar(bout::boundary::BoundaryRegionFCI* region, BoutReal value)
      : bndry(region), real_value(value), value_type(ValueType::REAL) {}
  BoundaryOpPar(bout::boundary::BoundaryRegionFCI* region)
      : bndry(region), real_value(0.), value_type(ValueType::REAL) {}
  BoundaryOpPar(bout::boundary::BoundaryRegionX* region,
                std::shared_ptr<FieldGenerator> value)
      : bndryX(region), gen_values(std::move(value)), value_type(ValueType::GEN) {}
  BoundaryOpPar(bout::boundary::BoundaryRegionX* region, Field3D* value)
      : bndryX(region), field_values(value), value_type(ValueType::FIELD) {}
  BoundaryOpPar(bout::boundary::BoundaryRegionX* region, BoutReal value)
      : bndryX(region), real_value(value) {}
  BoundaryOpPar(bout::boundary::BoundaryRegionX* region) : bndryX(region) {}
  BoundaryOpPar(bout::boundary::BoundaryRegionY* region,
                std::shared_ptr<FieldGenerator> value)
      : bndryY(region), gen_values(std::move(value)), value_type(ValueType::GEN) {}
  BoundaryOpPar(bout::boundary::BoundaryRegionY* region, Field3D* value)
      : bndryY(region), field_values(value), value_type(ValueType::FIELD) {}
  BoundaryOpPar(bout::boundary::BoundaryRegionY* region, BoutReal value)
      : bndryY(region), real_value(value) {}
  BoundaryOpPar(bout::boundary::BoundaryRegionY* region) : bndryY(region) {}
  BoundaryOpPar(BoundaryOpPar* region, std::shared_ptr<FieldGenerator> value)
      : bndry(region->bndry), gen_values(std::move(value)), value_type(ValueType::GEN) {}
  BoundaryOpPar(BoundaryOpPar* region, Field3D* value)
      : bndry(region->bndry), field_values(value), value_type(ValueType::FIELD) {}
  BoundaryOpPar(BoundaryOpPar* region, BoutReal value)
      : bndry(region->bndry), real_value(value) {}
  BoundaryOpPar(BoundaryOpPar* region) : bndry(region->bndry) {}
  ~BoundaryOpPar() override = default;

  // Note: All methods must implement clone, except for modifiers (see below)
  virtual BoundaryOpPar* clone(bout::boundary::BoundaryRegionFCI* region,
                               const std::list<std::string>& args) = 0;
  virtual BoundaryOpPar* clone(bout::boundary::BoundaryRegionFCI* region, Field3D* f) = 0;
  virtual BoundaryOpPar*
  clone(bout::boundary::BoundaryRegionFCI* region, const std::list<std::string>& args,
        const std::map<std::string, std::string>& UNUSED(keywords)) {
    return clone(region, args);
  }
  virtual BoundaryOpPar* clone(bout::boundary::BoundaryRegionX* region,
                               const std::list<std::string>& args) = 0;
  virtual BoundaryOpPar* clone(bout::boundary::BoundaryRegionX* region, Field3D* f) = 0;
  virtual BoundaryOpPar*
  clone(bout::boundary::BoundaryRegionX* region, const std::list<std::string>& args,
        const std::map<std::string, std::string>& UNUSED(keywords)) {
    return clone(region, args);
  }
  virtual BoundaryOpPar* clone(bout::boundary::BoundaryRegionY* region,
                               const std::list<std::string>& args) = 0;
  virtual BoundaryOpPar* clone(bout::boundary::BoundaryRegionY* region, Field3D* f) = 0;
  virtual BoundaryOpPar*
  clone(bout::boundary::BoundaryRegionY* region, const std::list<std::string>& args,
        const std::map<std::string, std::string>& UNUSED(keywords)) {
    return clone(region, args);
  }
  virtual BoundaryOpPar* clone(BoundaryOpPar* region,
                               const std::list<std::string>& args) = 0;
  virtual BoundaryOpPar* clone(BoundaryOpPar* region, Field3D* f) = 0;
  virtual BoundaryOpPar*
  clone(BoundaryOpPar* region, const std::list<std::string>& args,
        const std::map<std::string, std::string>& UNUSED(keywords)) {
    return clone(region, args);
  }

private:
  bout::boundary::BoundaryRegionFCI* bndry{nullptr};
  bout::boundary::BoundaryRegionX* bndryX{nullptr};
  bout::boundary::BoundaryRegionY* bndryY{nullptr};

  /// Possible ways to get boundary values
  std::shared_ptr<FieldGenerator> gen_values;
  Field3D* field_values{nullptr};
  BoutReal real_value{0.};

  /// Where to take boundary values from - the generator, field or BoutReal
  enum class ValueType { GEN, FIELD, REAL };
  const ValueType value_type{ValueType::REAL};

  BoutReal getValue(const bout::boundary::BoundaryRegionFCI::Point& bndry, BoutReal t);
  BoutReal getValue(const bout::boundary::BoundaryRegionX::Iterator& bndry, BoutReal t);
  BoutReal getValue(const bout::boundary::BoundaryRegionY::Iterator& bndry, BoutReal t);

  template <class T, bool isNeumann>
  friend class BoundaryOpParTemp;
};

template <class T, bool isNeumann = false>
class BoundaryOpParTemp : public BoundaryOpPar {
public:
  using BoundaryOpPar::BoundaryOpPar;

  using BoundaryOpPar::clone;

  BoundaryOpPar* clone(BoundaryOpPar* region,
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
  BoundaryOpPar* clone(bout::boundary::BoundaryRegionFCI* region,
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

  BoundaryOpPar* clone(bout::boundary::BoundaryRegionFCI* region, Field3D* f) override {
    return new T(region, f);
  }
  BoundaryOpPar* clone(bout::boundary::BoundaryRegionX* region,
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

  BoundaryOpPar* clone(bout::boundary::BoundaryRegionX* region, Field3D* f) override {
    return new T(region, f);
  }
  BoundaryOpPar* clone(bout::boundary::BoundaryRegionY* region,
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

  BoundaryOpPar* clone(bout::boundary::BoundaryRegionY* region, Field3D* f) override {
    return new T(region, f);
  }
  BoundaryOpPar* clone(BoundaryOpPar* region, Field3D* f) override {
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
    if (bndry != nullptr) {
      f.ynext(bndry->dir()).allocate(); // Ensure unique before modifying
      auto dy = f.getCoordinates()->dy();
      for (auto pnt : *bndry) {
        BoutReal value = getValue(pnt, t);
        if (isNeumann) {
          value *= dy[pnt.ind()];
        }
        static_cast<T*>(this)->apply_stencil(f, pnt, value);
      }
    }
    if (bndryX != nullptr) {
      f.allocate();
      auto dy = f.getCoordinates()->dx();
      for (auto pnt : *bndryX) {
        BoutReal value = getValue(pnt, t);
        if (isNeumann) {
          value *= dy[pnt.ind()];
        }
        static_cast<T*>(this)->apply_stencil(f, pnt, value);
      }
    }
    if (bndryY != nullptr) {
      f.allocate();
      auto dy = f.getCoordinates()->dy();
      for (auto pnt : *bndryY) {
        BoutReal value = getValue(pnt, t);
        if (isNeumann) {
          value *= dy[pnt.ind()];
        }
        static_cast<T*>(this)->apply_stencil(f, pnt, value);
      }
    }
  }
};

//////////////////////////////////////////////////
// Implementations

class BoundaryOpPar_dirichlet_o1 : public BoundaryOpParTemp<BoundaryOpPar_dirichlet_o1> {
public:
  using BoundaryOpParTemp::BoundaryOpParTemp;
  template <bout::boundary::BoundaryIterator Iter>
  static void apply_stencil(Field3D& f, Iter& pnt, BoutReal value) {
    dirichlet_o1(pnt, f, value);
  }
};

class BoundaryOpPar_dirichlet_o2 : public BoundaryOpParTemp<BoundaryOpPar_dirichlet_o2> {
public:
  using BoundaryOpParTemp::BoundaryOpParTemp;
  template <bout::boundary::BoundaryIterator Iter>
  static void apply_stencil(Field3D& f, Iter& pnt, BoutReal value) {
    dirichlet_o2(pnt, f, value);
  }
};

class BoundaryOpPar_dirichlet_o3 : public BoundaryOpParTemp<BoundaryOpPar_dirichlet_o3> {
public:
  using BoundaryOpParTemp::BoundaryOpParTemp;
  template <bout::boundary::BoundaryIterator Iter>
  static void apply_stencil(Field3D& f, Iter& pnt, BoutReal value) {
    dirichlet_o3(pnt, f, value);
  }
};

class BoundaryOpPar_neumann_o1
    : public BoundaryOpParTemp<BoundaryOpPar_neumann_o1, true> {
public:
  using BoundaryOpParTemp::BoundaryOpParTemp;
  template <bout::boundary::BoundaryIterator Iter>
  static void apply_stencil(Field3D& f, Iter& pnt, BoutReal value) {
    neumann_o1(pnt, f, value);
  }
};

class BoundaryOpPar_neumann_o2
    : public BoundaryOpParTemp<BoundaryOpPar_neumann_o2, true> {
public:
  using BoundaryOpParTemp::BoundaryOpParTemp;
  template <bout::boundary::BoundaryIterator Iter>
  static void apply_stencil(Field3D& f, Iter& pnt, BoutReal value) {
    neumann_o2(pnt, f, value);
  }
};

class BoundaryOpPar_neumann_o3
    : public BoundaryOpParTemp<BoundaryOpPar_neumann_o3, true> {
public:
  using BoundaryOpParTemp::BoundaryOpParTemp;
  template <bout::boundary::BoundaryIterator Iter>
  static void apply_stencil(Field3D& f, Iter& pnt, BoutReal value) {
    neumann_o3(pnt, f, value);
  }
};

#endif // BOUT_PAR_BNDRY_OP_H
