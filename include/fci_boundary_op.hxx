#ifndef __FCI_BNDRY_OP_H__
#define __FCI_BNDRY_OP_H__

#include <fci_derivs.hxx>
#include <boundary_op.hxx>

//////////////////////////////////////////////////
// Base class

class BoundaryOpFCI : public BoundaryOpBase {
public:
  BoundaryOpFCI();
  BoundaryOpFCI(BoundaryRegionFCI *region, FieldGenerator* value) :
    bndry(region),
    gen_values(value),
    value_type(GEN) {}
  BoundaryOpFCI(BoundaryRegionFCI *region, Field3D* value) :
    bndry(region),
    field_values(value),
    value_type(FIELD) {}
  BoundaryOpFCI(BoundaryRegionFCI *region, BoutReal value) :
    bndry(region),
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

  BoundaryRegionFCI *bndry;

protected:

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
  BoundaryOpFCI_dirichlet(BoundaryRegionFCI *region, FieldGenerator* value) :
    BoundaryOpFCI(region, value) {}
  BoundaryOpFCI_dirichlet(BoundaryRegionFCI *region, Field3D* value) :
    BoundaryOpFCI(region, value) {}
  BoundaryOpFCI_dirichlet(BoundaryRegionFCI *region, BoutReal value) :
    BoundaryOpFCI(region, value) {}

  void apply(Field3D &f) {return apply(f, 0);}
  void apply(Field3D &f, BoutReal t);

};

class BoundaryOpFCI_neumann : public BoundaryOpFCI {
public:
  BoundaryOpFCI_neumann();
  BoundaryOpFCI_neumann(BoundaryRegionFCI *region, FieldGenerator* value) :
    BoundaryOpFCI(region, value) {}
  BoundaryOpFCI_neumann(BoundaryRegionFCI *region, Field3D* value) :
    BoundaryOpFCI(region, value) {}
  BoundaryOpFCI_neumann(BoundaryRegionFCI *region, BoutReal value) :
    BoundaryOpFCI(region, value) {}

  void apply(Field3D &f) {return apply(f, 0);}
  void apply(Field3D &f, BoutReal t);

};

#endif // __FCI_BNDRY_OP_H__
