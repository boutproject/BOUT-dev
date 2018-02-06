#include <globals.hxx>
#include <output.hxx>
#include <utils.hxx>

#include "loadmetric.hxx"

void LoadMetric(BoutReal Lnorm, BoutReal Bnorm) {
  // Load metric coefficients from the mesh
  Field2D Rxy, Bpxy, Btxy, hthe, sinty;
  GRID_LOAD5(Rxy, Bpxy, Btxy, hthe, sinty); // Load metrics
  
  // Checking for dpsi and qinty used in BOUT grids
  Field2D dx;
  if(!mesh->get(dx,   "dpsi")) {
    output << "\tUsing dpsi as the x grid spacing\n";
    mesh->coordinates()->dx = dx; // Only use dpsi if found
  }else {
    // dx will have been read already from the grid
    output << "\tUsing dx as the x grid spacing\n";
  }
  //Field2D qinty;
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
  mesh->coordinates()->dx /= SQ(Lnorm)*Bnorm;
  
  Bpxy /= Bnorm;
  Btxy /= Bnorm;
  mesh->coordinates()->Bxy  /= Bnorm;
  
  // Calculate metric components
  string ptstr;
  Options::getRoot()->get("mesh:paralleltransform", ptstr, "identity");
  // Convert to lower case for comparison
  ptstr = lowercase(ptstr);
  if(ptstr == "shifted") {
    sinty = 0.0;  // I disappears from metric
  }
  
  BoutReal sbp = 1.0; // Sign of Bp
  if(min(Bpxy, true) < 0.0)
    sbp = -1.0;
  
  mesh->coordinates()->g11 = pow(Rxy*Bpxy,2);
  mesh->coordinates()->g22 = 1.0 / pow(hthe,2);
  mesh->coordinates()->g33 = pow(sinty,2)*mesh->coordinates()->g11 + pow(mesh->coordinates()->Bxy,2)/mesh->coordinates()->g11;
  mesh->coordinates()->g12 = 0.0;
  mesh->coordinates()->g13 = -sinty*mesh->coordinates()->g11;
  mesh->coordinates()->g23 = -sbp*Btxy/(hthe*Bpxy*Rxy);
  
  mesh->coordinates()->J = hthe / Bpxy;
  
  mesh->coordinates()->g_11 = 1.0/mesh->coordinates()->g11 + pow(sinty*Rxy,2);
  mesh->coordinates()->g_22 = pow(mesh->coordinates()->Bxy*hthe/Bpxy,2);
  mesh->coordinates()->g_33 = Rxy*Rxy;
  mesh->coordinates()->g_12 = sbp*Btxy*hthe*sinty*Rxy/Bpxy;
  mesh->coordinates()->g_13 = sinty*Rxy*Rxy;
  mesh->coordinates()->g_23 = sbp*Btxy*hthe*Rxy/Bpxy;
  
  mesh->coordinates()->geometry();
}

