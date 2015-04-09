#include "bout/constants.hxx"
#include "bout/mesh.hxx"
#include "field_factory.hxx"
#include "globals.hxx"
#include "output.hxx"
#include "parallel_boundary_op.hxx"

BoutReal BoundaryOpFCI::getValue(int x, int y, int z, BoutReal t) {

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

// void BoundaryOpFCI::apply_ddt(Field3D &f)

//////////////////////////////////////////
// Dirichlet boundary

void BoundaryOpFCI_dirichlet::apply(Field3D &f, BoutReal t) {

  BoundaryRegionFCI* bndry_fci = static_cast<BoundaryRegionFCI*>(bndry);

  Field3D& f_next = f.ynext(bndry_fci->dir);

  Coordinates& coord = *(mesh->coordinates());

  // Loop over grid points If point is in boundary, then fill in
  // f_next such that the field would be VALUE on the boundary
  for (bndry_fci->first(); !bndry_fci->isDone(); bndry_fci->next()) {
    // temp variables for convenience
    int x = bndry_fci->x; int y = bndry_fci->y; int z = bndry_fci->z;

    // Generate the boundary value
    BoutReal value = getValue(x, y, z, t);

    // Scale the field and normalise to the desired value
    BoutReal y_prime = bndry_fci->length;
    BoutReal f2 = (f(x,y,z) - value) * (coord.dy(x, y) - y_prime) / y_prime;

    f_next(x, y+bndry_fci->dir, z) = value - f2;
  }

}

//////////////////////////////////////////
// Neumann boundary

void BoundaryOpFCI_neumann::apply(Field3D &f, BoutReal t) {

  BoundaryRegionFCI* bndry_fci = static_cast<BoundaryRegionFCI*>(bndry);

  Field3D& f_next = f.ynext(bndry_fci->dir);

  // If point is in boundary, then fill in f_next such that the derivative
  // would be VALUE on the boundary
  for (bndry_fci->first(); !bndry_fci->isDone(); bndry_fci->next()) {
    // temp variables for convience
    int x = bndry_fci->x; int y = bndry_fci->y; int z = bndry_fci->z;

    // Generate the boundary value
    BoutReal value = getValue(x, y, z, t);

    f_next(x, y+bndry_fci->dir, z) = 2.*value + f(x, y, z);
  }

}
