#ifndef __PAR_BNDRY_OP_H__
#define __PAR_BNDRY_OP_H__

#include "boundary_op.hxx"
#include "bout_types.hxx"
#include "parallel_boundary_region.hxx"
#include "utils.hxx"

//////////////////////////////////////////////////
// Base class

class BoundaryOpPar : public BoundaryOpBase {
public:
  BoundaryOpPar();
  BoundaryOpPar(BoundaryRegionPar *region, FieldGenerator* value) :
    bndry(region),
    gen_values(value),
    value_type(GEN) {}
  BoundaryOpPar(BoundaryRegionPar *region, Field3D* value) :
    bndry(region),
    field_values(value),
    value_type(FIELD) {}
  BoundaryOpPar(BoundaryRegionPar *region, BoutReal value) :
    bndry(region),
    real_value(value),
    value_type(REAL) {}

  void apply(Field2D &f)
  {
    throw BoutException("Can't apply parallel boundary conditions to Field2D!");
  }
  void apply(Field2D &f, BoutReal t)
  {
    throw BoutException("Can't apply parallel boundary conditions to Field2D!");
  }
  void apply(Field3D &f) {}
  void apply(Field3D &f, BoutReal t) {}

  // Apply to time derivative
  // Unlikely to be used?
  void apply_ddt(Field3D &f) {};

  BoundaryRegionPar *bndry;

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

class BoundaryOpPar_dirichlet : public BoundaryOpPar {
public:
  BoundaryOpPar_dirichlet();
  BoundaryOpPar_dirichlet(BoundaryRegionPar *region, FieldGenerator* value) :
    BoundaryOpPar(region, value) {}
  BoundaryOpPar_dirichlet(BoundaryRegionPar *region, Field3D* value) :
    BoundaryOpPar(region, value) {}
  BoundaryOpPar_dirichlet(BoundaryRegionPar *region, BoutReal value) :
    BoundaryOpPar(region, value) {}

  void apply(Field3D &f) {return apply(f, 0);}
  void apply(Field3D &f, BoutReal t);

};

class BoundaryOpPar_neumann : public BoundaryOpPar {
public:
  BoundaryOpPar_neumann();
  BoundaryOpPar_neumann(BoundaryRegionPar *region, FieldGenerator* value) :
    BoundaryOpPar(region, value) {}
  BoundaryOpPar_neumann(BoundaryRegionPar *region, Field3D* value) :
    BoundaryOpPar(region, value) {}
  BoundaryOpPar_neumann(BoundaryRegionPar *region, BoutReal value) :
    BoundaryOpPar(region, value) {}

  void apply(Field3D &f) {return apply(f, 0);}
  void apply(Field3D &f, BoutReal t);

};

#endif // __PAR_BNDRY_OP_H__
