#include "bout/constants.hxx"
#include "bout/mesh.hxx"
#include "field_factory.hxx"
#include "globals.hxx"
#include "output.hxx"
#include "parallel_boundary_op.hxx"

BoutReal BoundaryOpPar::getValue(int x, int y, int z, BoutReal t) {

  BoutReal xnorm;
  BoutReal ynorm;
  BoutReal znorm;

  BoutReal value;

  switch (value_type) {
  case GEN:
    // This works but doesn't quite do the right thing... should
    // generate value on the boundary, but that gives wrong
    // answer. This instead generates the value at the gridpoint
    xnorm = mesh->GlobalX(x);
    ynorm = mesh->GlobalY(y);
    znorm = ((BoutReal)(z))/(mesh->ngz-1);
    return gen_values->generate(xnorm, TWOPI*ynorm, TWOPI*znorm, t);
  case FIELD:
    value = (*field_values)(x,y,z);
    return value;
  case REAL:
    return real_value;
  }

}

BoutReal BoundaryOpPar::getValue(const BoundaryRegionPar &bndry, BoutReal t) {

  BoutReal xnorm;
  BoutReal ynorm;
  BoutReal znorm;

  BoutReal value;

  switch (value_type) {
  case GEN:
    // Need to use GlobalX, except with BoutReal as argument...
    xnorm = mesh->GlobalX(bndry.s_x);
    ynorm = mesh->GlobalY(bndry.s_y);
    znorm = bndry.s_z/(mesh->ngz-1);
    return gen_values->generate(xnorm, TWOPI*ynorm, TWOPI*znorm, t);
  case FIELD:
    value = (*field_values)(bndry.x,bndry.y,bndry.z);
    return value;
  case REAL:
    return real_value;
  }

}

//////////////////////////////////////////
// Dirichlet boundary

BoundaryOpPar* BoundaryOpPar_dirichlet::clone(BoundaryRegionPar *region, const list<string> &args) {
  if(!args.empty()) {
    try {
      real_value = stringToReal(args.front());
      return new BoundaryOpPar_dirichlet(region, real_value);
    } catch (BoutException e) {
      FieldGenerator* newgen = 0;
      // First argument should be an expression
      newgen = FieldFactory::get()->parse(args.front());
      return new BoundaryOpPar_dirichlet(region, newgen);
    }
  }
  return new BoundaryOpPar_dirichlet(region);
}

BoundaryOpPar* BoundaryOpPar_dirichlet::clone(BoundaryRegionPar *region, Field3D *f) {
  return new BoundaryOpPar_dirichlet(region, f);
}

void BoundaryOpPar_dirichlet::apply(Field3D &f, BoutReal t) {

  Field3D& f_next = f.ynext(bndry->dir);

  Coordinates& coord = *(mesh->coordinates());

  // Loop over grid points If point is in boundary, then fill in
  // f_next such that the field would be VALUE on the boundary
  for (bndry->first(); !bndry->isDone(); bndry->next()) {
    // temp variables for convenience
    int x = bndry->x; int y = bndry->y; int z = bndry->z;

    // Generate the boundary value
    BoutReal value = getValue(*bndry, t);

    // Scale the field and normalise to the desired value
    BoutReal y_prime = bndry->length;
    BoutReal f2 = (f(x,y,z) - value) * (coord.dy(x, y) - y_prime) / y_prime;

    f_next(x, y+bndry->dir, z) = value - f2;
  }

}

//////////////////////////////////////////
// Neumann boundary

BoundaryOpPar* BoundaryOpPar_neumann::clone(BoundaryRegionPar *region, const list<string> &args) {
  if(!args.empty()) {
    try {
      real_value = stringToReal(args.front());
      return new BoundaryOpPar_neumann(region, real_value);
    } catch (BoutException e) {
      FieldGenerator* newgen = 0;
      // First argument should be an expression
      newgen = FieldFactory::get()->parse(args.front());
      return new BoundaryOpPar_neumann(region, newgen);
    }
  }
  return new BoundaryOpPar_neumann(region);
}

BoundaryOpPar* BoundaryOpPar_neumann::clone(BoundaryRegionPar *region, Field3D *f) {
  return new BoundaryOpPar_neumann(region, f);
}

void BoundaryOpPar_neumann::apply(Field3D &f, BoutReal t) {

  Field3D& f_next = f.ynext(bndry->dir);

  // If point is in boundary, then fill in f_next such that the derivative
  // would be VALUE on the boundary
  for (bndry->first(); !bndry->isDone(); bndry->next()) {
    // temp variables for convience
    int x = bndry->x; int y = bndry->y; int z = bndry->z;

    // Generate the boundary value
    BoutReal value = getValue(x, y, z, t);

    f_next(x, y+bndry->dir, z) = 2.*value + f(x, y, z);
  }

}
