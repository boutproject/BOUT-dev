#include <bout/tokamak_coordinates.hxx>

#include "loadmetric.hxx"

void LoadMetric(BoutReal Lnorm, BoutReal Bnorm) {
  // Load metric coefficients from the mesh
  Field2D Rxy, Bpxy, Btxy, hthe, sinty;

  auto mesh = bout::globals::mesh;
  auto coords = mesh->getCoordinates();

  GRID_LOAD5(Rxy, Bpxy, Btxy, hthe, sinty); // Load metrics

  // Checking for dpsi and qinty used in BOUT grids
  Field2D dx;
  if (!mesh->get(dx, "dpsi")) {
    output << "\tUsing dpsi as the x grid spacing\n";
    coords->setDx(dx); // Only use dpsi if found
  } else {
    // dx will have been read already from the grid
    output << "\tUsing dx as the x grid spacing\n";
  }
  Field2D qinty;

  const auto tokamak_coordinates_factory = TokamakCoordinatesFactory(*mesh);
  coords = tokamak_coordinates_factory.make_tokamak_coordinates();

  tokamak_coordinates_factory.normalise(Lnorm, Bnorm);

  // Calculate metric components
  std::string ptstr;
  Options::getRoot()->get("mesh:paralleltransform", ptstr, "identity");
  // Convert to lower case for comparison
  ptstr = lowercase(ptstr);
  if (ptstr == "shifted") {
    sinty = 0.0; // I disappears from metric
  }
}
