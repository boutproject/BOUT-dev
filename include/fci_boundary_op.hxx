#ifndef __FCI_BNDRY_OP_H__
#define __FCI_BNDRY_OP_H__

#include <fci_derivs.hxx>
#include <boundary_op.hxx>


class BoundaryFCI_dirichlet : public BoundaryOp {
  /// Possible ways to get boundary values
  FieldGenerator* gen_values;
  Field3D* field_values;
  BoutReal real_value;

  /// Where to take boundary values from - the generator, field or BoutReal
  enum ValueType { GEN, FIELD, REAL };
  const ValueType value_type;

  const FCIMap& fcimap;
  Field3D& f_next;

  BoutReal getValue(int x, int y, int z, BoutReal t);

  // Private default constructor
  BoundaryFCI_dirichlet();
public:
  BoundaryFCI_dirichlet(BoundaryRegionFCI* region, const FCIMap& fcimap, Field3D& f_next, FieldGenerator* value) :
    BoundaryOp(region),
    fcimap(fcimap),
    f_next(f_next),
    gen_values(value),
    value_type(GEN) {}
  BoundaryFCI_dirichlet(BoundaryRegion* region, const FCIMap& fcimap, Field3D& f_next, Field3D* value) :
    BoundaryOp(region),
    fcimap(fcimap),
    f_next(f_next),
    field_values(value),
    value_type(FIELD) {}
  BoundaryFCI_dirichlet(BoundaryRegion* region, const FCIMap& fcimap, Field3D& f_next, BoutReal value) :
    BoundaryOp(region),
    fcimap(fcimap),
    f_next(f_next),
    real_value(value),
    value_type(REAL) {}
  BoundaryOp* clone(BoundaryRegion* region, const list<string> &args);

  void apply(Field2D &f);
  void apply(Field2D &f, BoutReal t);
  void apply(Field3D &f);
  void apply(Field3D &f, BoutReal t);

  void apply_ddt(Field3D &f);

  BoundaryRegionFCI* bndry;
};

#endif // __FCI_BNDRY_OP_H__
