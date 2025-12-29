#include "bout/sys/generator_context.hxx"
#include "bout/boundary_region.hxx"
#include "bout/constants.hxx"
#include "bout/mesh.hxx"

namespace bout {
namespace generator {

Context::Context(int ix, int iy, int iz, CELL_LOC loc, Mesh* msh, BoutReal t)
    : ix_(ix), jy_(iy), kz_(iz), localmesh(msh) {

  parameters["x"] = (loc == CELL_XLOW) ? 0.5 * (msh->GlobalX(ix) + msh->GlobalX(ix - 1))
                                       : msh->GlobalX(ix);

  parameters["y"] = (loc == CELL_YLOW) ? PI * (msh->GlobalY(iy) + msh->GlobalY(iy - 1))
                                       : TWOPI * msh->GlobalY(iy);

  parameters["z"] = (loc == CELL_ZLOW)
                        ? TWOPI * (iz - 0.5) / static_cast<BoutReal>(msh->LocalNz)
                        : TWOPI * iz / static_cast<BoutReal>(msh->LocalNz);

  parameters["t"] = t;
}

Context::Context(const BoundaryRegion* bndry, int iz, CELL_LOC loc, BoutReal t, Mesh* msh)
    : // Add one to X index if boundary is in -x direction, so that XLOW is on the boundary
      ix_((bndry->bx < 0) ? bndry->x + 1 : bndry->x),
      jy_((bndry->by < 0) ? bndry->y + 1 : bndry->y), kz_(iz), localmesh(msh) {

  parameters["x"] = ((loc == CELL_XLOW) || (bndry->bx != 0))
                        ? 0.5 * (msh->GlobalX(ix_) + msh->GlobalX(ix_ - 1))
                        : msh->GlobalX(ix_);

  parameters["y"] = ((loc == CELL_YLOW) || (bndry->by != 0))
                        ? PI * (msh->GlobalY(jy_) + msh->GlobalY(jy_ - 1))
                        : TWOPI * msh->GlobalY(jy_);

  parameters["z"] = (loc == CELL_ZLOW)
                        ? TWOPI * (iz - 0.5) / static_cast<BoutReal>(msh->LocalNz)
                        : TWOPI * iz / static_cast<BoutReal>(msh->LocalNz);

  parameters["t"] = t;
}

} // namespace generator
} // namespace bout
