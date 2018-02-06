/*
 * Test LaplaceXZ solver
 *
 * Matrix assembly information:
 *    -mat_view ::ascii_info
 *
 * Useful diagnostics for SuperLU_dist solver
 * (pctype=lu, factor_package=superlu_dist)
 * -mat_superlu_dist_statprint
 */
#include <bout.hxx>

#include <bout/invert/laplacexy.hxx>
#include <field_factory.hxx>

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  ///////////////////////////////////////
  bool calc_metric;
  OPTION(Options::getRoot(), calc_metric, true);
  if(calc_metric) {
    // Read metric tensor
    Field2D Rxy, Btxy, Bpxy, B0, hthe, I;
    mesh->get(Rxy,  "Rxy");  // m
    mesh->get(Btxy, "Btxy"); // T
    mesh->get(Bpxy, "Bpxy"); // T
    mesh->get(B0,   "Bxy");  // T
    mesh->get(hthe, "hthe"); // m
    mesh->get(I,    "sinty");// m^-2 T^-1

    // Checking for dpsi used in BOUT grids
    Field2D dx;
    if(!mesh->get(dx,   "dpsi")) {
      output << "\tUsing dpsi as the x grid spacing\n";
      mesh->coordinates()->dx = dx; // Only use dpsi if found
    }else {
      // dx will have been read already from the grid
      output << "\tUsing dx as the x grid spacing\n";
    }

    Coordinates *coord = mesh->coordinates();

    // Calculate metric components
    string ptstr;
    Options::getRoot()->getSection("mesh")->get("paralleltransform", ptstr, "identity");
    // Convert to lower case for comparison
    ptstr = lowercase(ptstr);
    if(ptstr == "shifted") {
      I = 0.0;  // I disappears from metric
    }
      
    BoutReal sbp = 1.0; // Sign of Bp
    if(min(Bpxy, true) < 0.0)
      sbp = -1.0;

    // Calculate metrics
    coord->g11 = SQ(Rxy*Bpxy);
    coord->g22 = 1.0 / SQ(hthe);
    coord->g33 = SQ(I)*coord->g11 + SQ(B0)/coord->g11;
    coord->g12 = 0.0;
    coord->g13 = -I*coord->g11;
    coord->g23 = -sbp*Btxy/(hthe*Bpxy*Rxy);
    
    coord->J = hthe / Bpxy;
    coord->Bxy = B0;
    
    coord->g_11 = 1.0/coord->g11 + SQ(I*Rxy);
    coord->g_22 = SQ(B0*hthe/Bpxy);
    coord->g_33 = Rxy*Rxy;
    coord->g_12 = sbp*Btxy*hthe*I*Rxy/Bpxy;
    coord->g_13 = I*Rxy*Rxy;
    coord->g_23 = Btxy*hthe*Rxy/Bpxy;
    
    coord->geometry();
  }
  ///////////////////////////////////////
 
  LaplaceXY *inv = new LaplaceXY(mesh);

  output.write("Setting coefficients\n");

  Field2D A = FieldFactory::get()->create2D("a1", Options::getRoot(), mesh);
  Field2D B = FieldFactory::get()->create2D("b1", Options::getRoot(), mesh);
  //Field2D A = 1.;
  //Field2D B = 0.;
  inv->setCoefs(A,B);

  output.write("First solve\n");

  Field2D rhs = FieldFactory::get()->create2D("rhs", Options::getRoot(), mesh);
  Field2D x = inv->solve(rhs, 0.0);
  mesh->communicate(x);
  Field2D check = Laplace_perpXY(A,x) + B*x - rhs;
  Field2D check_laplaceperp = A*Laplace_perp(x) + B*x - rhs;

  SAVE_ONCE4(rhs, x, check, check_laplaceperp);

  output.write("Second solve\n");

  A = FieldFactory::get()->create2D("a2", Options::getRoot(), mesh);
  B = FieldFactory::get()->create2D("b2", Options::getRoot(), mesh);
  //A = 2.;
  //B = 0.1;
  inv->setCoefs(A,B);

  Field2D rhs2 = FieldFactory::get()->create2D("rhs2", Options::getRoot(), mesh);
  Field2D x2 = inv->solve(rhs2, 0.0);
  mesh->communicate(x2);
  Field2D check2 = Laplace_perpXY(A,x2) + B*x2 - rhs2;
  Field2D check2_laplaceperp = A*Laplace_perp(x2) + B*x2 - rhs2;
  SAVE_ONCE4(rhs2, x2, check2, check2_laplaceperp);

  dump.write();

  delete inv;

  BoutFinalise();
  return 0;
}
