#include <bout/tokamak_coordinates.hxx>

#include "loadmetric.hxx"

void LoadMetric(BoutReal Lnorm, BoutReal Bnorm) {

  auto mesh = bout::globals::mesh;
  auto coords = mesh->getCoordinates();

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

  // Calculate metric components
  std::string ptstr;
  Options::getRoot()->get("mesh:paralleltransform", ptstr, "identity");
  // Convert to lower case for comparison
  ptstr = lowercase(ptstr);
  BoutReal shearFactor = 1.0;
  if (ptstr == "shifted") {
    shearFactor = 0.0; // I disappears from metric
  }

  auto tokamak_options = TokamakOptions(*mesh);
  set_tokamak_coordinates_on_mesh(tokamak_options, *mesh, Lnorm, Bnorm, shearFactor);
}
