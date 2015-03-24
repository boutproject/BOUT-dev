#include <bout/physicsmodel.hxx>
#include <derivs.hxx>

class FieldAlign : public PhysicsModel {
protected:
  int init(bool restarting) {
    mesh->get(vx, "vx");
    mesh->get(vy, "vy");
    mesh->get(vz, "vz");
    mesh->get(G, "G");
    solver->add(f, "f");
    return 0;
  }
  
  int rhs(BoutReal t) {
    	mesh->communicate(f);
	f.applyBoundary(t);

	ddt(f) = 
		vx / G * (mesh->g11*DDX(f) + mesh->g12*DDY(f) + mesh->g13*DDZ(f)) +
		vy / G * (mesh->g12*DDX(f) + mesh->g22*DDY(f) + mesh->g23*DDZ(f)) +    // Upwinding with second-order central differencing 
		vz / G * (mesh->g13*DDX(f) + mesh->g23*DDY(f) + mesh->g33*DDZ(f));      // (unstable without additional dissipation)
		- SQ(SQ(mesh->dx))*D4DX4(f) - SQ(SQ(mesh->dy))*D4DY4(f) - SQ(SQ(mesh->dz))*D4DY4(f);   // Numerical dissipation terms
    
    
    return 0;
  }
  
private:
  Field3D f;
  Field2D vx, vy, vz, G;
};

BOUTMAIN(FieldAlign);
