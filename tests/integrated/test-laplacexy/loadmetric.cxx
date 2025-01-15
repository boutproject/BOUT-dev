#include <bout/tokamak_coordinates.hxx>

#include "loadmetric.hxx"

void LoadMetric(BoutReal Lnorm, BoutReal Bnorm) {

  bool noshear;

  auto mesh = bout::globals::mesh;
  auto coords = mesh->getCoordinates();

  // Calculate metric components
  std::string ptstr;
  Options::getRoot()->get("mesh:paralleltransform", ptstr, "identity");
  // Convert to lower case for comparison
  ptstr = lowercase(ptstr);
  if (ptstr == "shifted") {
    noshear = true;
  }

  auto tokamak_options = TokamakOptions(*mesh);
  set_tokamak_coordinates_on_mesh(tokamak_options, *mesh, true, Lnorm, Bnorm);
}
