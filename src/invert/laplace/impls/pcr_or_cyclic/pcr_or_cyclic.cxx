#include "pcr_or_cyclic.hxx"

#include "bout/mesh.hxx"
#include "output.hxx"

#include "../cyclic/cyclic_laplace.hxx"
#include "../pcr/pcr.hxx"

LaplacePCRorCyclic::LaplacePCRorCyclic(Options* opt, const CELL_LOC loc, Mesh* mesh_in) {
  const Mesh& localmesh = (mesh_in == nullptr) ? *bout::globals::mesh : *mesh_in;

  // Number of X procs must be a power of 2
  const bool x_procs_pow2 = is_pow2(localmesh.getNXPE());
  // Number of x points must be a power of 2
  const bool x_points_pow2 = is_pow2(localmesh.GlobalNxNoBoundaries);

  if (x_procs_pow2 and x_points_pow2) {
    output_info.write("Using LaplacePCR as internal Laplacian solver\n");
    pcr_or_cyclic = std::make_unique<LaplacePCR>(opt, loc, mesh_in);
  } else {
    output_info.write("Using LaplaceCyclic as internal Laplacian solver\n");
    output_warn.write(
        "WARNING: Unable to use 'pcr' Laplacian inversion solver, falling back to\n"
        "         'cyclic'! This has reduced performance, especially at large core\n"
        "         counts. 'pcr' requires the number of processors in x to be a\n"
        "         power-of-two (currently this is {}) AND the number of points in x\n"
        "         (excluding boundaries) to be a power-of-two (currently this is {})\n",
        localmesh.getNXPE(), localmesh.GlobalNxNoBoundaries);
    pcr_or_cyclic = std::make_unique<LaplaceCyclic>(opt, loc, mesh_in);
  }
}
