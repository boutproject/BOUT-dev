#ifndef __FCI_BNDRY_OP_H__
#define __FCI_BNDRY_OP_H__

#include <fci_derivs.hxx>
#include <boundary_op.hxx>

//////////////////////////////////////////////////
// Base class

class BoundaryOpFCI : public BoundaryOp {
public:
  BoundaryOpFCI();
  BoundaryOpFCI(const FCIMap& fcimap, FieldGenerator* value) :
    BoundaryOp(fcimap.boundary),
    fcimap(fcimap),
    gen_values(value),
    value_type(GEN) {}
  BoundaryOpFCI(const FCIMap& fcimap, Field3D* value) :
    BoundaryOp(fcimap.boundary),
    fcimap(fcimap),
    field_values(value),
    value_type(FIELD) {}
  BoundaryOpFCI(const FCIMap& fcimap, BoutReal value) :
    BoundaryOp(fcimap.boundary),
    fcimap(fcimap),
    real_value(value),
    value_type(REAL) {}

  void apply(Field2D &f)
  {
      throw BoutException("Can't apply FCI boundary conditions to Field2D!");
  }
  void apply(Field2D &f, BoutReal t)
  {
      throw BoutException("Can't apply FCI boundary conditions to Field2D!");
  }
  void apply(Field3D &f) {}
  void apply(Field3D &f, BoutReal t) {}

  // Apply to time derivative
  // Unlikely to be used?
  void apply_ddt(Field3D &f) {};

protected:
  const FCIMap& fcimap;

  /// Possible ways to get boundary values
  FieldGenerator* gen_values;
  Field3D* field_values;
  BoutReal real_value;

  /// Where to take boundary values from - the generator, field or BoutReal
  enum ValueType { GEN, FIELD, REAL };
  const ValueType value_type;

  BoutReal getValue(int x, int y, int z, BoutReal t);

};

//////////////////////////////////////////////////
// Implementations

class BoundaryOpFCI_dirichlet : public BoundaryOpFCI {
public:
  BoundaryOpFCI_dirichlet();
  BoundaryOpFCI_dirichlet(const FCIMap& fcimap, FieldGenerator* value) :
    BoundaryOpFCI(fcimap, value) {}
  BoundaryOpFCI_dirichlet(const FCIMap& fcimap, Field3D* value) :
    BoundaryOpFCI(fcimap, value) {}
  BoundaryOpFCI_dirichlet(const FCIMap& fcimap, BoutReal value) :
    BoundaryOpFCI(fcimap, value) {}

  void apply(Field3D &f) {return apply(f, 0);}
  void apply(Field3D &f, BoutReal t);

};

class BoundaryOpFCI_neumann : public BoundaryOpFCI {
public:
  BoundaryOpFCI_neumann();
  BoundaryOpFCI_neumann(const FCIMap& fcimap, FieldGenerator* value) :
    BoundaryOpFCI(fcimap, value) {}
  BoundaryOpFCI_neumann(const FCIMap& fcimap, Field3D* value) :
    BoundaryOpFCI(fcimap, value) {}
  BoundaryOpFCI_neumann(const FCIMap& fcimap, BoutReal value) :
    BoundaryOpFCI(fcimap, value) {}

  void apply(Field3D &f) {return apply(f, 0);}
  void apply(Field3D &f, BoutReal t);

};

#endif // __FCI_BNDRY_OP_H__
