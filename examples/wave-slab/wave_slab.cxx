/*
  Simple wave test in a sheared slab domain.
  
  Uses the same field-aligned Clebsch coordinate system as
  most BOUT++ tokamak simulations. See the coordinates manual
  for details.
  
  Note: Here the only components of the coordinate system which 
  are tested are g_22 (for Grad_par), and the twist shift angle.
  
 */

#include <bout/physicsmodel.hxx>


class WaveTest : public PhysicsModel {
public:
  int init(bool restarting) {
    
    Field2D Rxy, Bpxy, Btxy, hthe, I;
    GRID_LOAD(Rxy);
    GRID_LOAD(Bpxy);
    GRID_LOAD(Btxy);
    GRID_LOAD(hthe);
    mesh->get(mesh->Bxy,  "Bxy");
    
    if(mesh->ShiftXderivs) {
      // No integrated shear in metric
      I = 0.0;
    }else
      mesh->get(I,    "sinty");
    
    mesh->g11 = (Rxy*Bpxy)^2;
    mesh->g22 = 1.0 / (hthe^2);
    mesh->g33 = (I^2)*mesh->g11 + (mesh->Bxy^2)/mesh->g11;
    mesh->g12 = 0.0;
    mesh->g13 = -I*mesh->g11;
    mesh->g23 = -Btxy/(hthe*Bpxy*Rxy);
    
    mesh->J = hthe / Bpxy;
    
    mesh->g_11 = 1.0/mesh->g11 + ((I*Rxy)^2);
    mesh->g_22 = (mesh->Bxy*hthe/Bpxy)^2;
    mesh->g_33 = Rxy*Rxy;
    mesh->g_12 = Btxy*hthe*I*Rxy/Bpxy;
    mesh->g_13 = I*Rxy*Rxy;
    mesh->g_23 = Btxy*hthe*Rxy/Bpxy;
    
    mesh->geometry();
    
    solver->add(f, "f");
    solver->add(g, "g");
    
    return 0;
  }
  int rhs(BoutReal time) {
    mesh->communicate(f,g);
    
    ddt(f) = Grad_par(g);
    ddt(g) = Grad_par(f);
    
    return 0;
  }
  
private:
  Field3D f, g;
};


BOUTMAIN(WaveTest);
