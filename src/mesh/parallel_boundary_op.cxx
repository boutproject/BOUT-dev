#include "parallel_boundary_op.hxx"
#include "field_factory.hxx"
#include "globals.hxx"
#include "output.hxx"
#include "bout/constants.hxx"
#include "bout/mesh.hxx"

BoutReal BoundaryOpPar::getValue(const BoundaryRegionPar& bndry, BoutReal t) {
  BoutReal value;

  switch (value_type) {
  case ValueType::GEN:
    return gen_values->generate(bout::generator::Context(
        bndry.s_x(), bndry.s_y(), bndry.s_z(), CELL_CENTRE, bndry.localmesh, t));
  case ValueType::FIELD:
    // FIXME: Interpolate to s_x, s_y, s_z...
    value = (*field_values)[bndry.ind()];
    return value;
  case ValueType::REAL:
    return real_value;
  default:
    throw BoutException("Invalid value_type encountered in BoundaryOpPar::getValue");
  }
}
