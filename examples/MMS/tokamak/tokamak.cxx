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
  int init(bool restarting) {
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
      - 10.*(SQ(SQ(mesh->dx))*D4DX4(advect) + SQ(SQ(mesh->dz))*D4DZ4(advect));
    
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
  
    // Checking for dpsi and qinty used in BOUT grids
    Field2D dx;
    if(!mesh->get(dx,   "dpsi")) {
      output << "\tUsing dpsi as the x grid spacing\n";
      mesh->dx = dx; // Only use dpsi if found
    }else {
      // dx will have been read already from the grid
      output << "\tUsing dx as the x grid spacing\n";
    }
    Field2D qinty;
    if(!mesh->get(qinty, "qinty")) {
      output << "\tUsing qinty as the Z shift\n";
      mesh->zShift = qinty;
    }else {
      // Keep zShift
      output << "\tUsing zShift as the Z shift\n";
    }

    Rxy      /= Lnorm;
    hthe     /= Lnorm;
    sinty    *= SQ(Lnorm)*Bnorm;
    mesh->dx /= SQ(Lnorm)*Bnorm;
  
    Bpxy /= Bnorm;
    Btxy /= Bnorm;
    mesh->Bxy  /= Bnorm;
  
    // Calculate metric components
    if(mesh->ShiftXderivs) {
      sinty = 0.0;  // I disappears from metric
    }
  
    BoutReal sbp = 1.0; // Sign of Bp
    if(min(Bpxy, true) < 0.0)
      sbp = -1.0;

    mesh->g11 = (Rxy*Bpxy)^2;
    mesh->g22 = 1.0 / (hthe^2);
    mesh->g33 = (sinty^2)*mesh->g11 + (mesh->Bxy^2)/mesh->g11;
    mesh->g12 = 0.0;
    mesh->g13 = -sinty*mesh->g11;
    mesh->g23 = -sbp*Btxy/(hthe*Bpxy*Rxy);
  
    mesh->J = hthe / Bpxy;
  
    mesh->g_11 = 1.0/mesh->g11 + ((sinty*Rxy)^2);
    mesh->g_22 = (mesh->Bxy*hthe/Bpxy)^2;
    mesh->g_33 = Rxy*Rxy;
    mesh->g_12 = sbp*Btxy*hthe*sinty*Rxy/Bpxy;
    mesh->g_13 = sinty*Rxy*Rxy;
    mesh->g_23 = sbp*Btxy*hthe*Rxy/Bpxy;
  
    mesh->geometry();
  }

private:
  Field3D drive;
  Field3D laplacepar, delp2, advect;
};

BOUTMAIN(TokamakMMS);
