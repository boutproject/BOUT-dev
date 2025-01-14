#include <bout/tokamak_coordinates.hxx>

#include "loadmetric.hxx"

void LoadMetric(BoutReal Lnorm, BoutReal Bnorm) {

  bool noshear;

  auto mesh = bout::globals::mesh;

  // Calculate metric components
  std::string ptstr;
  Options::getRoot()->get("mesh:paralleltransform", ptstr, "identity");
  // Convert to lower case for comparison
  ptstr = lowercase(ptstr);
  if (ptstr == "shifted") {
    noshear = true;
  }

  auto tokamak_coordinates = TokamakCoordinates(*mesh);
  set_tokamak_coordinates_on_mesh(tokamak_coordinates, *mesh, true, Lnorm, Bnorm);
}
