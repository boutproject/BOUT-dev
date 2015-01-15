#include <fci_boundary.hxx>

void BoundaryRegionFCI::add_point(const int x, const int y, const int z) {
  bndry_points.push_back({x, y, z});
}

void BoundaryRegionFCI::first() {
  bndry_position = bndry_points.begin();
  if (!isDone()) {
	x = bndry_position->x;
	y = bndry_position->y;
	z = bndry_position->z;
  }
}

void BoundaryRegionFCI::next() {
  ++bndry_position;
  if (!isDone()) {
	x = bndry_position->x;
	y = bndry_position->y;
	z = bndry_position->z;
  }
}

bool BoundaryRegionFCI::isDone() {
  return (bndry_position == bndry_points.end());
}

BoutReal BoundaryFCI_dirichlet::getValue(int x, int y, int z) {

  switch (value_type) {
  case GEN:
	// This works but doesn't quite do the right thing... should
	// generate value on the boundary, but that gives wrong
	// answer. This instead generates the value at the gridpoint
	BoutReal xnorm = mesh->GlobalX(x);
	BoutReal ynorm = mesh->GlobalY(y);
	BoutReal znorm = ((BoutReal)(z))/(mesh->ngz-1);
	return gen_values->generate(xnorm, TWOPI*ynorm, TWOPI*znorm, t);
  case FIELD:
    return field_values[x][y][z];
  case REAL:
    return real_value;
  }   

}

void BoundaryFCI_dirichlet::apply(Field3D &f) {

  // Loop over grid points If point is in boundary, then fill in
  // f_next such that the field would be VALUE on the boundary
  for (bndry->first(); !bndry->isDone(); bndry->next()) {
    // temp variables for convience
	int x = bndry->x; int y = bndry->y; int	z = bndry->z;

	// Generate the boundary value
    getValue(x, y, z)

	// Scale the field and normalise to the desired value
	BoutReal y_prime = fcimap.y_prime[x][y][z];
	BoutReal f2 = (f[x][y][z] - value) * (mesh.dy(x, y) - y_prime) / y_prime;

	f_next[x][y+fcimap.dir][z] = value - f2;
  }

}

void BoundaryFCI::neumannBC(Field3D &f, Field3D &f_next, const FCIMap &fcimap) {
  // If point is in boundary, then fill in f_next such that the derivative
  // would be VALUE on the boundary
  for (bndry->first(); !bndry->isDone(); bndry->next()) {
    // temp variables for convience
	int x = bndry->x; int y = bndry->y; int	z = bndry->z;

	f_next[x][y+fcimap.dir][z] = f[x][y][z];
  }
}

