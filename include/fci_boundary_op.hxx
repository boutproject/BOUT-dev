#ifndef __FCI_BNDRY_OP_H__
#define __FCI_BNDRY_OP_H__

#include <fci_derivs.hxx>
#include <boundary_op.hxx>


class BoundaryFCI_dirichlet : public BoundaryOp {

  const FCIMap& fcimap;

  /// Possible ways to get boundary values
  FieldGenerator* gen_values;
  Field3D* field_values;
  BoutReal real_value;

  /// Where to take boundary values from - the generator, field or BoutReal
  enum ValueType { GEN, FIELD, REAL };
  const ValueType value_type;

  BoutReal getValue(int x, int y, int z, BoutReal t);

  // Private default constructor
  BoundaryFCI_dirichlet();
public:
  BoundaryFCI_dirichlet(const FCIMap& fcimap, FieldGenerator* value) :
    BoundaryOp(fcimap.boundary),
    fcimap(fcimap),
    gen_values(value),
    value_type(GEN) {}
  BoundaryFCI_dirichlet(const FCIMap& fcimap, Field3D* value) :
    BoundaryOp(fcimap.boundary),
    fcimap(fcimap),
    field_values(value),
    value_type(FIELD) {}
  BoundaryFCI_dirichlet(const FCIMap& fcimap, BoutReal value) :
    BoundaryOp(fcimap.boundary),
    fcimap(fcimap),
    real_value(value),
    value_type(REAL) {}
  BoundaryOp* clone(BoundaryRegion* region, const list<string> &args);

  void apply(Field2D &f);
  void apply(Field2D &f, BoutReal t);
  void apply(Field3D &f);
  void apply(Field3D &f, BoutReal t);

  void apply_ddt(Field3D &f);
};

#endif // __FCI_BNDRY_OP_H__
