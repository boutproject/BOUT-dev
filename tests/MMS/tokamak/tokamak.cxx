/*
 * MMS test in tokamak geometry
 *
 * Independently tests the bracket operator, perpendicular diffusion,
 * and parallel diffusion operators.
 * 
 */

#include <bout/physicsmodel.hxx>
#include <field_factory.hxx>
#include <derivs.hxx>

class TokamakMMS : public PhysicsModel {
public:
  int init(bool UNUSED(restarting)) {
    solver->add(laplacepar, "laplacepar");
    solver->add(delp2, "delp2");
    solver->add(advect, "advect");
    
    // Load the metric tensor
    LoadMetric(1.0, 1.0);

    return 0;
  }
  int rhs(BoutReal time) {
    mesh->communicate(advect, delp2, laplacepar);
    
    drive = FieldFactory::get()->create3D("drive:solution", Options::getRoot(), mesh, CELL_CENTRE, time);
    
    // Test bracket advection operator
    ddt(advect) = -1e-3*bracket(drive, advect, BRACKET_ARAKAWA)
      - 10.*(SQ(SQ(mesh->getCoordinates()->dx))*D4DX4(advect) + SQ(SQ(mesh->getCoordinates()->dz))*D4DZ4(advect));
    
    // Test perpendicular diffusion operator
    ddt(delp2) = 1e-5*Delp2(delp2);
    
    // Test parallel diffusion operator
    ddt(laplacepar) = Laplace_par(laplacepar);
    
    return 0;
  }
  void LoadMetric(BoutReal Lnorm, BoutReal Bnorm) {
    // Load metric coefficients from the mesh
    Field2D Rxy, Bpxy, Btxy, hthe, sinty;
    GRID_LOAD5(Rxy, Bpxy, Btxy, hthe, sinty); // Load metrics

    Coordinates *coords = mesh->getCoordinates();
  
    // Checking for dpsi used in BOUT grids
    Field2D dx;
    if(!mesh->get(dx,   "dpsi")) {
      output << "\tUsing dpsi as the x grid spacing\n";
      coords->dx = dx; // Only use dpsi if found
    }else {
      // dx will have been read already from the grid
      output << "\tUsing dx as the x grid spacing\n";
    }

    Rxy      /= Lnorm;
    hthe     /= Lnorm;
    sinty    *= SQ(Lnorm)*Bnorm;
    coords->dx /= SQ(Lnorm)*Bnorm;
  
    Bpxy /= Bnorm;
    Btxy /= Bnorm;
    coords->Bxy /= Bnorm;
  
    // Calculate metric components
    bool ShiftXderivs;
    Options::getRoot()->get("shiftXderivs", ShiftXderivs, false); // Read global flag
    if(ShiftXderivs) {
     sinty = 0.0;  // I disappears from metric
    }
  
    BoutReal sbp = 1.0; // Sign of Bp
    if(min(Bpxy, true) < 0.0)
      sbp = -1.0;

    coords->g11 = SQ(Rxy*Bpxy);
    coords->g22 = 1.0 / SQ(hthe);
    coords->g33 = SQ(sinty)*coords->g11 + SQ(coords->Bxy)/coords->g11;
    coords->g12 = 0.0;
    coords->g13 = -sinty*coords->g11;
    coords->g23 = -sbp*Btxy/(hthe*Bpxy*Rxy);
  
    coords->J = hthe / Bpxy;
  
    coords->g_11 = 1.0/coords->g11 + SQ(sinty*Rxy);
    coords->g_22 = SQ(coords->Bxy*hthe/Bpxy);
    coords->g_33 = Rxy*Rxy;
    coords->g_12 = sbp*Btxy*hthe*sinty*Rxy/Bpxy;
    coords->g_13 = sinty*Rxy*Rxy;
    coords->g_23 = sbp*Btxy*hthe*Rxy/Bpxy;
  
    coords->geometry();
  }

private:
  Field3D drive;
  Field3D laplacepar, delp2, advect;
};

BOUTMAIN(TokamakMMS);
