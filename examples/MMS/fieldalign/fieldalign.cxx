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

	// df/dt = df/dpsi + df/dtheta + df/dphi
// 	output<<" "<<mesh->g12[0][1]<<" "<<mesh->g13[0][1]<<" "<<mesh->g23[0][1]<<" "<<mesh->g11[0][1]<<" "<<mesh->g22[0][1]<<" "<<mesh->g33[0][1]<<"\n";
// 	output<<" "<<vx[3][3][3]<<" "<<vy[3][3][3]<<" "<<vz[3][3][3]<<"\n";
// 	output<<DDX(f)[5][3][3]<<"\n";
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
