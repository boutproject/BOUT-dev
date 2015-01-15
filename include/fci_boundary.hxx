#ifndef __FCI_BNDRY_H__
#define __FCI_BNDRY_H__

#include <fci_derivs.hxx>
#include <boundary_op.hxx>
#include <boundary_region.hxx>

/**
 * Boundary region for FCI. This contains a vector of points that are
 * inside the boundary.
 *
 */
class BoundaryRegionFCI : public BoundaryRegion {

  struct Indices {
    int x;
    int y;
    int z;
  };

  typedef std::vector<Indices> IndicesVec;
  typedef IndicesVec::iterator IndicesIter;
  
  /// Vector of points in the boundary
  IndicesVec bndry_points;
  /// Current position in the boundary points
  IndicesIter bndry_position;
  
public:
  BoundaryRegionFCI(const string &name, BndryLoc loc) : label(name), location(loc) {}
  
  /// Add a point to the boundary
  void add_point(const int x, const int y, const int z);

  void first();
  void next();
  bool isDone();
  
  /// Index of the point in the boundary
  int z;
  
};

class BoundaryFCI_dirichlet : public BoundaryOp {
  /// Possible ways to get boundary values
  FieldGenerator* gen_values;
  Field3D* field_values;
  BoutReal real_value;

  /// Where to take boundary values from - the generator, field or BoutReal
  enum ValueType { GEN, FIELD, REAL };
  const ValueType value_type;

  FCIMap* fcimap;
  Field3D* f_next;

  BoutReal getValue(int x, int y, int z);

  // Private default constructor
  BoundaryFCI_dirichlet();
public:
  BoundaryFCI_dirichlet(BoundaryRegion* region, FCIMap* fcimap, Field3D* f_next, FieldGenerator* value) :
    bndry(region),
    fcimap(fcimap),
    f_next(f_next),
    gen_values(value),
    value_type(GEN) {}
  BoundaryFCI_dirichlet(BoundaryRegion* region, FCIMap* fcimap, Field3D* f_next, Field3D* value) :
    bndry(region),
    fcimap(fcimap),
    f_next(f_next),
    field_values(value),
    value_type(FIELD) {}
  BoundaryFCI_dirichlet(BoundaryRegion* region, FCIMap* fcimap, Field3D* f_next, BoutReal value) :
    bndry(region),
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
};

#endif //  __FCI_BNDRY_H__
