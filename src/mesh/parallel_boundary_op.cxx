#include "bout/constants.hxx"
#include "bout/mesh.hxx"
#include "field_factory.hxx"
#include "globals.hxx"
#include "output.hxx"
#include "parallel_boundary_op.hxx"

BoutReal BoundaryOpPar::getValue(int x, int y, int z, BoutReal t) {

  Mesh* mesh = bndry->localmesh;

  BoutReal xnorm;
  BoutReal ynorm;
  BoutReal znorm;

  BoutReal value;

  switch (value_type) {
  case ValueType::GEN:
    // This works but doesn't quite do the right thing... should
    // generate value on the boundary, but that gives wrong
    // answer. This instead generates the value at the gridpoint
    xnorm = mesh->GlobalX(x);
    ynorm = mesh->GlobalY(y);
    znorm = static_cast<BoutReal>(z) / (mesh->LocalNz);
    return gen_values->generate(xnorm, TWOPI*ynorm, TWOPI*znorm, t);
  case ValueType::FIELD:
    value = (*field_values)(x,y,z);
    return value;
  case ValueType::REAL:
    return real_value;
  default:
    throw BoutException("Invalid value_type encountered in BoundaryOpPar::getValue");
  }

}

BoutReal BoundaryOpPar::getValue(const BoundaryRegionPar &bndry, BoutReal t) {

  Mesh* mesh = bndry.localmesh;

  BoutReal xnorm;
  BoutReal ynorm;
  BoutReal znorm;

  BoutReal value;

  switch (value_type) {
  case ValueType::GEN:
    // Need to use GlobalX, except with BoutReal as argument...
    xnorm = mesh->GlobalX(bndry.s_x);
    ynorm = mesh->GlobalY(bndry.s_y);
    znorm = bndry.s_z/(mesh->LocalNz);
    return gen_values->generate(xnorm, TWOPI*ynorm, TWOPI*znorm, t);
  case ValueType::FIELD:
    // FIXME: Interpolate to s_x, s_y, s_z...
    value = (*field_values)(bndry.x,bndry.y,bndry.z);
    return value;
  case ValueType::REAL:
    return real_value;
  default:
    throw BoutException("Invalid value_type encountered in BoundaryOpPar::getValue");
  }

}

//////////////////////////////////////////
// Dirichlet boundary

BoundaryOpPar* BoundaryOpPar_dirichlet::clone(BoundaryRegionPar *region, const std::list<std::string> &args) {
  if(!args.empty()) {
    try {
      real_value = stringToReal(args.front());
      return new BoundaryOpPar_dirichlet(region, real_value);
    } catch (BoutException& e) {
      std::shared_ptr<FieldGenerator> newgen = nullptr;
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

  Coordinates& coord = *(f.getCoordinates());

  // Loop over grid points If point is in boundary, then fill in
  // f_next such that the field would be VALUE on the boundary
  for (bndry->first(); !bndry->isDone(); bndry->next()) {
    // temp variables for convenience
    int x = bndry->x; int y = bndry->y; int z = bndry->z;

    // Generate the boundary value
    BoutReal value = getValue(*bndry, t);

    // Scale the field and normalise to the desired value
    BoutReal y_prime = bndry->length;
    BoutReal f2 = (f(x, y, z) - value) * (coord.dy(x, y, z) - y_prime) / y_prime;

    f_next(x, y+bndry->dir, z) = value - f2;
  }
}

//////////////////////////////////////////
// Dirichlet boundary - Third order

BoundaryOpPar* BoundaryOpPar_dirichlet_O3::clone(BoundaryRegionPar *region, const std::list<std::string> &args) {
  if(!args.empty()) {
    try {
      real_value = stringToReal(args.front());
      return new BoundaryOpPar_dirichlet_O3(region, real_value);
    } catch (BoutException& e) {
      std::shared_ptr<FieldGenerator> newgen = nullptr;
      // First argument should be an expression
      newgen = FieldFactory::get()->parse(args.front());
      return new BoundaryOpPar_dirichlet_O3(region, newgen);
    }
  }
  return new BoundaryOpPar_dirichlet_O3(region);
}

BoundaryOpPar* BoundaryOpPar_dirichlet_O3::clone(BoundaryRegionPar *region, Field3D *f) {
  return new BoundaryOpPar_dirichlet_O3(region, f);
}

void BoundaryOpPar_dirichlet_O3::apply(Field3D &f, BoutReal t) {

  Field3D& f_next = f.ynext(bndry->dir);
  Field3D& f_prev = f.ynext(-bndry->dir);

  Coordinates& coord = *(f.getCoordinates());

  // Loop over grid points If point is in boundary, then fill in
  // f_next such that the field would be VALUE on the boundary
  for (bndry->first(); !bndry->isDone(); bndry->next()) {
    // temp variables for convenience
    int x = bndry->x; int y = bndry->y; int z = bndry->z;

    // Generate the boundary value
    BoutReal fb = getValue(*bndry, t);
    BoutReal f1 = f_prev(x, y-bndry->dir, z);
    BoutReal f2 = f(x,y,z);
    BoutReal l1 = coord.dy(x, y, z);
    BoutReal l2 = bndry->length;
    BoutReal l3 = coord.dy(x, y, z) - l2;

    BoutReal denom = (l1*l1*l2 + l1*l2*l2);
    BoutReal term1 = (l2*l2*l3 + l2*l3*l3);
    BoutReal term2 = l1*(l1+l2+l3)*(l2+l3);
    BoutReal term3 = l3*((l1+l2)*l3 + (l1+l2)*(l1+l2));

    f_next(x, y+bndry->dir, z) = (term1*f1 + term2*fb - term3*f2)/denom;
  }
}

//////////////////////////////////////////
// Dirichlet with interpolation

BoundaryOpPar* BoundaryOpPar_dirichlet_interp::clone(BoundaryRegionPar *region, const std::list<std::string> &args) {
  if(!args.empty()) {
    try {
      real_value = stringToReal(args.front());
      return new BoundaryOpPar_dirichlet_interp(region, real_value);
    } catch (BoutException& e) {
      std::shared_ptr<FieldGenerator> newgen = nullptr;
      // First argument should be an expression
      newgen = FieldFactory::get()->parse(args.front());
      return new BoundaryOpPar_dirichlet_interp(region, newgen);
    }
  }
  return new BoundaryOpPar_dirichlet_interp(region);
}

BoundaryOpPar* BoundaryOpPar_dirichlet_interp::clone(BoundaryRegionPar *region, Field3D *f) {
  return new BoundaryOpPar_dirichlet_interp(region, f);
}

void BoundaryOpPar_dirichlet_interp::apply(Field3D &f, BoutReal t) {

  Field3D& f_next = f.ynext(bndry->dir);
  Field3D& f_prev = f.ynext(-bndry->dir);

  Coordinates& coord = *(f.getCoordinates());

  // Loop over grid points If point is in boundary, then fill in
  // f_next such that the field would be VALUE on the boundary
  for (bndry->first(); !bndry->isDone(); bndry->next()) {
    // temp variables for convenience
    int x = bndry->x; int y = bndry->y; int z = bndry->z;

    // Generate the boundary value
    BoutReal fs = getValue(*bndry, t);

    // Scale the field and normalise to the desired value
    BoutReal dy = coord.dy(x, y, z);
    BoutReal s = bndry->length*dy;

    f_next(x, y+bndry->dir, z) = f_prev(x, y-bndry->dir, z)*(1.-(2.*s/(dy+s)))
      + 2.*f(x, y, z)*((s-dy)/s)
      + fs*(dy/s - (2./s + 1.));
  }

}

//////////////////////////////////////////
// Neumann boundary

BoundaryOpPar* BoundaryOpPar_neumann::clone(BoundaryRegionPar *region, const std::list<std::string> &args) {
  if(!args.empty()) {
    try {
      real_value = stringToReal(args.front());
      return new BoundaryOpPar_neumann(region, real_value);
    } catch (BoutException& e) {
      std::shared_ptr<FieldGenerator> newgen = nullptr;
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
  TRACE("BoundaryOpPar_neumann::apply");
  
  Field3D& f_next = f.ynext(bndry->dir);
  f_next.allocate(); // Ensure unique before modifying
  
  Coordinates& coord = *(f.getCoordinates());

  // If point is in boundary, then fill in f_next such that the derivative
  // would be VALUE on the boundary
  for (bndry->first(); !bndry->isDone(); bndry->next()) {
    // temp variables for convience
    int x = bndry->x; int y = bndry->y; int z = bndry->z;

    // Generate the boundary value
    BoutReal value = getValue(x, y, z, t);
    BoutReal dy = coord.dy(x, y, z);

    f_next(x, y+bndry->dir, z) = f(x, y, z) + bndry->dir*value*dy;
  }

}
