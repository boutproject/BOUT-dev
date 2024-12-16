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
  auto coords = tokamak_coordinates.make_coordinates(noshear, Lnorm, Bnorm);
}
