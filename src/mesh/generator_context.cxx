#include "bout/sys/generator_context.hxx"
#include "boundary_region.hxx"
#include "bout/constants.hxx"
#include "bout/mesh.hxx"

Context::Context(int ix, int iy, int iz, CELL_LOC loc, Mesh* msh, BoutReal t) : localmesh(msh) {
  
  parameters["x"] = (loc == CELL_XLOW) ?
    0.5 * (msh->GlobalX(ix) + msh->GlobalX(ix - 1)) :
    msh->GlobalX(ix);
  
  parameters["y"] = (loc == CELL_YLOW) ?
    PI * (msh->GlobalY(iy) + msh->GlobalY(iy - 1)) :
    TWOPI * msh->GlobalY(iy);

  parameters["z"] = (loc == CELL_ZLOW) ?
    TWOPI * (iz - 0.5) / static_cast<BoutReal>(msh->LocalNz) :
    TWOPI * iz / static_cast<BoutReal>(msh->LocalNz);

  parameters["t"] = t;
}

Context::Context(const BoundaryRegion* bndry, int iz, CELL_LOC loc, BoutReal t, Mesh * msh) : localmesh(msh) {
  
  // Add one to X index if boundary is in -x direction, so that XLOW is on the boundary
  int ix = (bndry->bx < 0) ? bndry->x + 1 : bndry->x;
  
  parameters["x"] = ((loc == CELL_XLOW) || (bndry->bx != 0)) ?
    0.5 * (msh->GlobalX(ix) + msh->GlobalX(ix - 1)) :
    msh->GlobalX(ix);

  int iy = (bndry->by < 0) ? bndry->y + 1 : bndry->y;

  parameters["y"] = ((loc == CELL_YLOW) || bndry->by) ?
    PI * (msh->GlobalY(iy) + msh->GlobalY(iy - 1)) :
    TWOPI * msh->GlobalY(iy);

  parameters["z"] = (loc == CELL_ZLOW) ?
    TWOPI * (iz - 0.5) / static_cast<BoutReal>(msh->LocalNz) :
    TWOPI * iz / static_cast<BoutReal>(msh->LocalNz);

  parameters["t"] = t;
}
