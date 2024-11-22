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

  auto tokamak_coordinates_factory = TokamakCoordinatesFactory(*mesh);
  tokamak_coordinates_factory.normalise(Lnorm, Bnorm);
  auto coords = tokamak_coordinates_factory.make_tokamak_coordinates(noshear);
}
