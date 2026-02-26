#include "bout/parallel_boundary_op.hxx"
#include "bout/bout_types.hxx"
#include "bout/constants.hxx"
#include "bout/field_factory.hxx"
#include "bout/globals.hxx"
#include "bout/mesh.hxx"
#include "bout/output.hxx"
#include "bout/parallel_boundary_region.hxx"

BoutReal BoundaryOpPar::getValue(const BoundaryRegionParIter& bndry, BoutReal t) {
  switch (value_type) {
  case ValueType::GEN:
    return gen_values->generate(bout::generator::Context(
        bndry.s_x(), bndry.s_y(), bndry.s_z(), CELL_CENTRE, bndry.localmesh, t));
  case ValueType::FIELD:
    // FIXME: Interpolate to s_x, s_y, s_z...
    return (*field_values)[bndry.ind()];
  case ValueType::REAL:
    return real_value;
  default:
    throw BoutException("Invalid value_type encountered in BoundaryOpPar::getValue");
  }
}
