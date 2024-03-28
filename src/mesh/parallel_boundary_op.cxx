#include "bout/parallel_boundary_op.hxx"
#include "bout/constants.hxx"
#include "bout/field_factory.hxx"
#include "bout/globals.hxx"
#include "bout/mesh.hxx"
#include "bout/output.hxx"

BoutReal BoundaryOpPar::getValue(const BoundaryRegionPar& bndry, BoutReal t) {

  Mesh* mesh = bndry.localmesh;

  BoutReal value;

  switch (value_type) {
  case ValueType::GEN:
    return gen_values->generate(
        bout::generator::Context(bndry.s_x, bndry.s_y, bndry.s_z, CELL_CENTRE, mesh, t));
  case ValueType::FIELD:
    // FIXME: Interpolate to s_x, s_y, s_z...
    value = (*field_values)(bndry.x, bndry.y, bndry.z);
    return value;
  case ValueType::REAL:
    return real_value;
  default:
    throw BoutException("Invalid value_type encountered in BoundaryOpPar::getValue");
  }
}

//////////////////////////////////////////
// Dirichlet boundary

void BoundaryOpPar_dirichlet::apply(Field3D& f, BoutReal t) {
  Field3D& f_next = f.ynext(bndry->dir);

  Coordinates const& coord = *(f.getCoordinates());

  // Loop over grid points If point is in boundary, then fill in
  // f_next such that the field would be VALUE on the boundary
  for (bndry->first(); !bndry->isDone(); bndry->next()) {
    // temp variables for convenience
    int const x = bndry->x;
    int const y = bndry->y;
    int const z = bndry->z;

    // Generate the boundary value
    BoutReal const value = getValue(*bndry, t);

    // Scale the field and normalise to the desired value
    BoutReal const y_prime = bndry->length;
    BoutReal const f2 = (f(x, y, z) - value) * (coord.dy()(x, y, z) - y_prime) / y_prime;

    f_next(x, y + bndry->dir, z) = value - f2;
  }
}

//////////////////////////////////////////
// Dirichlet boundary - Third order

void BoundaryOpPar_dirichlet_O3::apply(Field3D& f, BoutReal t) {

  Field3D& f_next = f.ynext(bndry->dir);
  Field3D& f_prev = f.ynext(-bndry->dir);

  Coordinates const& coord = *(f.getCoordinates());

  // Loop over grid points If point is in boundary, then fill in
  // f_next such that the field would be VALUE on the boundary
  for (bndry->first(); !bndry->isDone(); bndry->next()) {
    // temp variables for convenience
    int const x = bndry->x;
    int const y = bndry->y;
    int const z = bndry->z;

    // Generate the boundary value
    BoutReal const fb = getValue(*bndry, t);
    BoutReal const f1 = f_prev(x, y - bndry->dir, z);
    BoutReal const f2 = f(x, y, z);
    BoutReal const l1 = coord.dy()(x, y, z);
    BoutReal const l2 = bndry->length;
    BoutReal const l3 = coord.dy()(x, y, z) - l2;

    BoutReal const denom = (l1 * l1 * l2 + l1 * l2 * l2);
    BoutReal const term1 = (l2 * l2 * l3 + l2 * l3 * l3);
    BoutReal const term2 = l1 * (l1 + l2 + l3) * (l2 + l3);
    BoutReal const term3 = l3 * ((l1 + l2) * l3 + (l1 + l2) * (l1 + l2));

    f_next(x, y + bndry->dir, z) = (term1 * f1 + term2 * fb - term3 * f2) / denom;
  }
}

//////////////////////////////////////////
// Dirichlet with interpolation

void BoundaryOpPar_dirichlet_interp::apply(Field3D& f, BoutReal t) {

  Field3D& f_next = f.ynext(bndry->dir);
  Field3D& f_prev = f.ynext(-bndry->dir);

  Coordinates const& coord = *(f.getCoordinates());

  // Loop over grid points If point is in boundary, then fill in
  // f_next such that the field would be VALUE on the boundary
  for (bndry->first(); !bndry->isDone(); bndry->next()) {
    // temp variables for convenience
    int const x = bndry->x;
    int const y = bndry->y;
    int const z = bndry->z;

    // Generate the boundary value
    BoutReal const fs = getValue(*bndry, t);

    // Scale the field and normalise to the desired value
    BoutReal const dy = coord.dy()(x, y, z);
    BoutReal const s = bndry->length * dy;

    f_next(x, y + bndry->dir, z) =
        f_prev(x, y - bndry->dir, z) * (1. - (2. * s / (dy + s)))
        + 2. * f(x, y, z) * ((s - dy) / s) + fs * (dy / s - (2. / s + 1.));
  }
}

//////////////////////////////////////////
// Neumann boundary

void BoundaryOpPar_neumann::apply(Field3D& f, BoutReal t) {
  TRACE("BoundaryOpPar_neumann::apply");

  Field3D& f_next = f.ynext(bndry->dir);
  f_next.allocate(); // Ensure unique before modifying

  Coordinates const& coord = *(f.getCoordinates());

  // If point is in boundary, then fill in f_next such that the derivative
  // would be VALUE on the boundary
  for (bndry->first(); !bndry->isDone(); bndry->next()) {
    // temp variables for convience
    int const x = bndry->x;
    int const y = bndry->y;
    int const z = bndry->z;

    // Generate the boundary value
    BoutReal const value = getValue(*bndry, t);
    BoutReal const dy = coord.dy()(x, y, z);

    f_next(x, y + bndry->dir, z) = f(x, y, z) + bndry->dir * value * dy;
  }
}
