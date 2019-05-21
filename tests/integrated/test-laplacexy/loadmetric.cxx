#include "bout/mesh.hxx"
#include "field2d.hxx"
#include "globals.hxx"
#include "output.hxx"
#include "utils.hxx"

#include "loadmetric.hxx"

void LoadMetric(BoutReal Lnorm, BoutReal Bnorm) {
  // Load metric coefficients from the mesh
  Field2D Rxy, Bpxy, Btxy, hthe, sinty;

  auto mesh = bout::globals::mesh;
  auto coords = mesh->getCoordinates();

  GRID_LOAD5(Rxy, Bpxy, Btxy, hthe, sinty); // Load metrics
  
  // Checking for dpsi and qinty used in BOUT grids
  Field2D dx;
  if(!mesh->get(dx,   "dpsi")) {
    output << "\tUsing dpsi as the x grid spacing\n";
    coords->dx = dx; // Only use dpsi if found
  }else {
    // dx will have been read already from the grid
    output << "\tUsing dx as the x grid spacing\n";
  }
  Field2D qinty;
  //if(!mesh->get(qinty, "qinty")) {
  //  output << "\tUsing qinty as the Z shift\n";
  //  mesh->getParallelTransform().zShift = qinty;
  //}else {
  //  // Keep zShift
  //  output << "\tUsing zShift as the Z shift\n";
  //}

  Rxy      /= Lnorm;
  hthe     /= Lnorm;
  sinty    *= SQ(Lnorm)*Bnorm;
  coords->dx /= SQ(Lnorm)*Bnorm;
  
  Bpxy /= Bnorm;
  Btxy /= Bnorm;
  coords->Bxy  /= Bnorm;
  
  // Calculate metric components
  std::string ptstr;
  Options::getRoot()->get("mesh:paralleltransform", ptstr, "identity");
  // Convert to lower case for comparison
  ptstr = lowercase(ptstr);
  if(ptstr == "shifted") {
    sinty = 0.0;  // I disappears from metric
  }
  
  BoutReal sbp = 1.0; // Sign of Bp
  if(min(Bpxy, true) < 0.0)
    sbp = -1.0;
  
  coords->g11 = pow(Rxy*Bpxy,2);
  coords->g22 = 1.0 / pow(hthe,2);
  coords->g33 = pow(sinty,2)*coords->g11 + pow(coords->Bxy,2)/coords->g11;
  coords->g12 = 0.0;
  coords->g13 = -sinty*coords->g11;
  coords->g23 = -sbp*Btxy/(hthe*Bpxy*Rxy);
  
  coords->J = hthe / Bpxy;
  
  coords->g_11 = 1.0/coords->g11 + pow(sinty*Rxy,2);
  coords->g_22 = pow(coords->Bxy*hthe/Bpxy,2);
  coords->g_33 = Rxy*Rxy;
  coords->g_12 = sbp*Btxy*hthe*sinty*Rxy/Bpxy;
  coords->g_13 = sinty*Rxy*Rxy;
  coords->g_23 = sbp*Btxy*hthe*Rxy/Bpxy;
  
  coords->geometry();
}

