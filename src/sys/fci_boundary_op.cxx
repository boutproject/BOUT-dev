#include <fci_boundary_op.hxx>
#include <output.hxx>

BoundaryOp* BoundaryFCI_dirichlet::clone(BoundaryRegion* region, const list<string> &args) {

  return this;
  // return new BoundaryFCI_dirichlet(region);

}

BoutReal BoundaryFCI_dirichlet::getValue(int x, int y, int z, BoutReal t) {

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
    value = (*field_values)[x][y][z];
    return value;
  case REAL:
    return real_value;
  }   

}

void BoundaryFCI_dirichlet::apply(Field2D &f) {
  // error
  output << "Can't apply FCI boundary conditions to Field2D!\n";
}
  
void BoundaryFCI_dirichlet::apply(Field2D &f, BoutReal t) {
  // error
  output << "Can't apply FCI boundary conditions to Field2D!\n";
}


void BoundaryFCI_dirichlet::apply(Field3D &f) {

  apply(f, 0);

}

void BoundaryFCI_dirichlet::apply(Field3D &f, BoutReal t) {

  Field3D* f_next;

  switch(bndry->location) {
  case BNDRY_FCI_FWD:
    f_next = f.yup();
    break;
  case BNDRY_FCI_BKWD:
    f_next = f.ydown();
    break;
  default:
    output << "Can't apply FCI boundary to non-FCI region\n";
    // throw exception
  }

  BoundaryRegionFCI* bndry_fci = static_cast<BoundaryRegionFCI*>(bndry);

  // Loop over grid points If point is in boundary, then fill in
  // f_next such that the field would be VALUE on the boundary
  for (bndry_fci->first(); !bndry_fci->isDone(); bndry_fci->next()) {
    // temp variables for convenience
	int x = bndry_fci->x; int y = bndry_fci->y; int	z = bndry_fci->z;

	// Generate the boundary value
    BoutReal value = getValue(x, y, z, t);

	// Scale the field and normalise to the desired value
	BoutReal y_prime = fcimap.y_prime[x][y][z];
	BoutReal f2 = (f[x][y][z] - value) * (mesh->dy(x, y) - y_prime) / y_prime;

	(*f_next)(x, y+fcimap.dir, z) = value - f2;
  }

}

// void BoundaryFCI_dirichlet::apply_ddt(Field2D &f) {
//   // error
//   output << "Can't apply FCI boundary conditions to Field2D!\n";
// }

void BoundaryFCI_dirichlet::apply_ddt(Field3D &f) {

  Field3D* dt = f.timeDeriv();

  apply(*dt, 0);

}

// void BoundaryFCI::neumannBC(Field3D &f, Field3D &f_next, const FCIMap &fcimap) {
//   // If point is in boundary, then fill in f_next such that the derivative
//   // would be VALUE on the boundary
//   for (bndry_fci->first(); !bndry_fci->isDone(); bndry_fci->next()) {
//     // temp variables for convience
// 	int x = bndry_fci->x; int y = bndry_fci->y; int	z = bndry_fci->z;

// 	f_next[x][y+fcimap->dir][z] = f[x][y][z];
//   }
// }

