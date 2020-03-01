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
  int init(bool UNUSED(restarting)) {
    auto *coords = mesh->getCoordinates();
    Field2D Rxy, Bpxy, Btxy, hthe, I;
    GRID_LOAD(Rxy);
    GRID_LOAD(Bpxy);
    GRID_LOAD(Btxy);
    GRID_LOAD(hthe);
    mesh->get(coords->Bxy,  "Bxy");
    int ShiftXderivs = 0;
    mesh->get(ShiftXderivs, "false");
    if(ShiftXderivs) {
      // No integrated shear in metric
      I = 0.0;
    }else
      mesh->get(I,    "sinty");
    
    coords->g11 = pow(Rxy*Bpxy,2.0);
    coords->g22 = 1.0 / pow(hthe,2.0);
    coords->g33 = pow(I,2.0)*coords->g11 + pow(coords->Bxy,2.0)/coords->g11;
    coords->g12 = 0.0;
    coords->g13 = -I*coords->g11;
    coords->g23 = -Btxy/(hthe*Bpxy*Rxy);
    
    coords->J = hthe / Bpxy;
    
    coords->g_11 = 1.0/coords->g11 + (pow(I*Rxy,2.0));
    coords->g_22 = pow(coords->Bxy*hthe/Bpxy,2.0);
    coords->g_33 = Rxy*Rxy;
    coords->g_12 = Btxy*hthe*I*Rxy/Bpxy;
    coords->g_13 = I*Rxy*Rxy;
    coords->g_23 = Btxy*hthe*Rxy/Bpxy;
    
    coords->geometry();
    
    solver->add(f, "f");
    solver->add(g, "g");
    
    return 0;
  }
  int rhs(BoutReal UNUSED(time)) {
    mesh->communicate(f,g);
    
    ddt(f) = Grad_par(g);
    ddt(g) = Grad_par(f);
    
    return 0;
  }
  
private:
  Field3D f, g;
};


BOUTMAIN(WaveTest);
