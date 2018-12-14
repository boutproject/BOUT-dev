#include "bout/sys/generator_context.hxx"
#include "bout/boundary_region.hxx"
#include "bout/constants.hxx"
#include "bout/mesh.hxx"

namespace bout {
namespace generator {

Context::Context(int ix, int iy, int iz, CELL_LOC loc, Mesh* msh, BoutReal t)
    : localmesh(msh) {

  parameters["x"] = (loc == CELL_XLOW) ? 0.5 * (msh->GlobalX(ix) + msh->GlobalX(ix - 1))
                                       : msh->GlobalX(ix);

  parameters["y"] = (loc == CELL_YLOW) ? PI * (msh->GlobalY(iy) + msh->GlobalY(iy - 1))
                                       : TWOPI * msh->GlobalY(iy);

  parameters["z"] = (loc == CELL_ZLOW) ? PI * (msh->GlobalZ(iz) + msh->GlobalZ(iz - 1))
                                       : TWOPI * msh->GlobalZ(iz);

  parameters["t"] = t;
}

Context::Context(const BoundaryRegion* bndry, int iz, CELL_LOC loc, BoutReal t, Mesh* msh)
    : localmesh(msh) {

  // Add one to X index if boundary is in -x direction, so that XLOW is on the boundary
  const int ix = (bndry->bx < 0) ? bndry->x + 1 : bndry->x;

  parameters["x"] = ((loc == CELL_XLOW) || (bndry->bx != 0))
                        ? 0.5 * (msh->GlobalX(ix) + msh->GlobalX(ix - 1))
                        : msh->GlobalX(ix);

  const int iy = (bndry->by < 0) ? bndry->y + 1 : bndry->y;

  parameters["y"] = ((loc == CELL_YLOW) || (bndry->by != 0))
                        ? PI * (msh->GlobalY(iy) + msh->GlobalY(iy - 1))
                        : TWOPI * msh->GlobalY(iy);

  parameters["z"] = (loc == CELL_ZLOW) ? PI * (msh->GlobalZ(iz) + msh->GlobalZ(iz - 1))
                                       : TWOPI * msh->GlobalZ(iz);

  parameters["t"] = t;
}

Context::Context(BoutReal x, BoutReal y, BoutReal z, Mesh* msh, BoutReal t)
    : localmesh(msh), parameters{{"x", x}, {"y", y}, {"z", z}, {"t", t}} {}

} // namespace generator
} // namespace bout
