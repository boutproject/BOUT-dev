#include <bout.hxx>

#include <bout/invert/laplacexy.hxx>
#include <derivs.hxx>
#include <field_factory.hxx>

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);
  
  ///////////////////////////////////////
  bool calc_metric;
  OPTION(Options::getRoot(), calc_metric, false);
  if(calc_metric) {
    // Read metric tensor
    Field2D Rxy, Btxy, Bpxy, B0, hthe, I;
    mesh->get(Rxy,  "Rxy");  // m
    mesh->get(Btxy, "Btxy"); // T
    mesh->get(Bpxy, "Bpxy"); // T
    mesh->get(B0,   "Bxy");  // T
    mesh->get(hthe, "hthe"); // m
    mesh->get(I,    "sinty");// m^-2 T^-1
    
    // Calculate metrics
    mesh->g11 = (Rxy*Bpxy)^2;
    mesh->g22 = 1.0 / (hthe^2);
    mesh->g33 = (I^2)*mesh->g11 + (B0^2)/mesh->g11;
    mesh->g12 = 0.0;
    mesh->g13 = -I*mesh->g11;
    mesh->g23 = -Btxy/(hthe*Bpxy*Rxy);
    
    mesh->J = hthe / Bpxy;
    mesh->Bxy = B0;
    
    mesh->g_11 = 1.0/mesh->g11 + ((I*Rxy)^2);
    mesh->g_22 = (B0*hthe/Bpxy)^2;
    mesh->g_33 = Rxy*Rxy;
    mesh->g_12 = Btxy*hthe*I*Rxy/Bpxy;
    mesh->g_13 = I*Rxy*Rxy;
    mesh->g_23 = Btxy*hthe*Rxy/Bpxy;
    
    mesh->geometry();
  }
  ///////////////////////////////////////
  
  // Read an analytic input
  Field2D input = FieldFactory::get()->create2D("input", Options::getRoot(), mesh);
  
  // Create a LaplaceXY solver
  LaplaceXY *laplacexy = new LaplaceXY(mesh);
  
  // Solve, using 0.0 as starting guess
  Field2D solved = laplacexy->solve(input, 0.0);
  
  // Need to communicate guard cells
  mesh->communicate(solved);
  
  // Now differentiate using Laplace_perp
  Field2D result = Laplace_perp(solved);

  // Write fields to output
  SAVE_ONCE3(input, solved, result);
  dump.write();
  
  BoutFinalise();
  return 0;
}

