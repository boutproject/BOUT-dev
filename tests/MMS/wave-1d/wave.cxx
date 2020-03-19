#include <bout/physicsmodel.hxx>
#include <initialprofiles.hxx>
#include <derivs.hxx>
#include <cmath>
#include "mathematica.h"
#include <bout/constants.hxx>
#include <unused.hxx>

BoutReal MS_f(BoutReal t, BoutReal  x, BoutReal  y, BoutReal  z);
BoutReal dxMS_f(BoutReal t, BoutReal  x, BoutReal  y, BoutReal  z);

BoutReal MS_g(BoutReal t, BoutReal  x, BoutReal  y, BoutReal  z);
BoutReal dxMS_g(BoutReal t, BoutReal  x, BoutReal  y, BoutReal  z);

BoutReal Lx, Ly, Lz; // Size of the domain

class Wave1D : public PhysicsModel {
private:
  Field3D f, g; // Evolving variables
  Field3D E_f, E_g; // Error vectors

  Coordinates *coord;
  
  const Field3D solution_f(BoutReal t);
  const Field3D source_f(BoutReal t);
  
  const Field3D solution_g(BoutReal t);
  const Field3D source_g(BoutReal t);
protected:
  int init(bool) override {
    // Coordinate system
    coord = mesh->getCoordinates();
    
    // Get the options
    Options &meshoptions = Options::root()["mesh"];
    
    Lx = meshoptions["Lx"].withDefault(1.0);
    Ly = meshoptions["Ly"].withDefault(1.0);
    
    // this assumes equidistant grid
    int nguard = mesh->xstart;
    coord->dx = Lx / (mesh->GlobalNx - 2 * nguard);
    coord->dy = Ly / (mesh->GlobalNy - 2 * nguard);

    SAVE_ONCE(Lx, Ly);
    
    //set mesh
    coord->g11 = 1.0;
    coord->g22 = 1.0;
    coord->g33 = 1.0;
    coord->g12 = 0.0;
    coord->g13 = 0.0;
    coord->g23 = 0.0;
    
    coord->g_11 = 1.0;
    coord->g_22 = 1.0;
    coord->g_33 = 1.0;
    coord->g_12 = 0.0;
    coord->g_13 = 0.0;
    coord->g_23 = 0.0;
    coord->geometry();

    g.setLocation(CELL_XLOW); // g staggered to the left of f
    
    // Dirichlet everywhere except inner x-boundary Neumann
    f.addBndryFunction(MS_f,BNDRY_ALL);
    f.addBndryFunction(dxMS_f,BNDRY_XIN);

    g.addBndryFunction(MS_g,BNDRY_ALL);
    g.addBndryFunction(dxMS_g,BNDRY_XIN);
    
    // Tell BOUT++ to solve f and g
    SOLVE_FOR(f, g);

    // Set initial condition to MS at t = 0.
    if (mesh->StaggerGrids) {
      for (int xi = mesh->xstart; xi < mesh->xend + 1; xi++) {
        for (int yj = mesh->ystart; yj < mesh->yend + 1; yj++){
	  for (int zk = 0; zk < mesh->LocalNz; zk++) {
	    f(xi, yj, zk) = MS_f(0.,mesh->GlobalX(xi),mesh->GlobalY(yj),coord->dz(xi,yj,zk)*zk);
	    g(xi, yj, zk) = MS_g(0.,0.5*(mesh->GlobalX(xi)+mesh->GlobalX(xi-1)),mesh->GlobalY(yj),coord->dz(xi,yj,zk)*zk);

            output.write("{:d}: {:e}\n", xi, g(xi, yj, zk));
          }
	}
      }
    } else {
      for (int xi = mesh->xstart; xi < mesh->xend + 1; xi++) {
        for (int yj = mesh->ystart; yj < mesh->yend + 1; yj++) {
          for (int zk = 0; zk < mesh->LocalNz; zk++) {
            f(xi, yj, zk) = MS_f(0., mesh->GlobalX(xi), mesh->GlobalY(yj), coord->dz(xi, yj, zk) * zk);
            g(xi, yj, zk) = MS_g(0., mesh->GlobalX(xi), mesh->GlobalY(yj), coord->dz(xi, yj, zk) * zk);
          }
        }
      }
    }
    E_f.allocate();
    E_g.allocate();
    SAVE_REPEAT(E_f, E_g);
    
    return 0;
  }

  int rhs(BoutReal t) override {
    mesh->communicate(f,g); // Communicate guard cells
    
    //update time-dependent boundary conditions
    f.applyBoundary(t);
    g.applyBoundary(t);

    // Central differencing
    ddt(f) = DDX(g, CELL_CENTRE);// + 20*SQ(coord->dx)*D2DX2(f);
    ddt(g) = DDX(f, CELL_XLOW);// + 20*SQ(coord->dx)*D2DX2(g);
    
    // Add MMS source term
    ddt(f) += source_f(t);
    ddt(g) += source_g(t);
    
    return 0;
  }

  // This called every output timestep
  int outputMonitor(BoutReal simtime, int UNUSED(iter), int UNUSED(NOUT)) override {

    Field3D Sf = solution_f(simtime);
    Field3D Sg = solution_g(simtime);
    
    //Calculate the error. norms are calculated in the post-processing
    E_f = 0.;
    E_g = 0.;
    for (int xi = mesh->xstart; xi < mesh->xend + 1; xi++){
      for (int yj = mesh->ystart ; yj < mesh->yend + 1; yj++){
        for (int zk = 0; zk < mesh->LocalNz ; zk++) {
          
          E_f(xi,yj,zk) = f(xi,yj,zk) - Sf(xi,yj,zk);
          E_g(xi,yj,zk) = g(xi,yj,zk) - Sg(xi,yj,zk);
          output_error.write("Error at {:d},{:d},{:d} = {:e}, {:e}",
                       xi, yj, zk, E_f(xi,yj,zk), E_g(xi,yj,zk));
        }
      }
    }
    return 0;
  }
};


/////////////////// SOLUTION FOR F //////////////////////////////

//Manufactured solution
BoutReal MS_f(BoutReal t, BoutReal  x, BoutReal  UNUSED(y), BoutReal  UNUSED(z)) {
  // Input is in normalised x,y,z location
  x *= Lx;         // X input [0,1]
  return (BoutReal)0.9 + 0.9*x + 0.2*Cos(10*t)*Sin(5.*Power(x,2));
}

//x-derivative of MS. For Neumann bnd cond
BoutReal dxMS_f(BoutReal t, BoutReal  x, BoutReal  UNUSED(y), BoutReal  UNUSED(z)) {
  x *= Lx;         // X input [0,1]
  return 0.9 + 2.*x*Cos(10*t)*Cos(5.*Power(x,2));
}

//Manufactured solution
const Field3D Wave1D::solution_f(BoutReal t) {
  Field3D S;
  S.allocate();
  
  int bx = (mesh->LocalNx - (mesh->xend - mesh->xstart + 1)) / 2;
  int by = (mesh->LocalNy - (mesh->yend - mesh->ystart + 1)) / 2;

  for (int xi = mesh->xstart - bx; xi < mesh->xend + bx + 1; xi++){
    for (int yj = mesh->ystart - by; yj < mesh->yend + by + 1; yj++){
      BoutReal x = mesh->GlobalX(xi);
      BoutReal y = mesh->GlobalY(yj);//GlobalY not fixed yet
      for (int zk = 0; zk < mesh->LocalNz; zk++) {
        BoutReal z = coord->dz(xi,yj,zk)*zk;
        S(xi, yj, zk) = MS_f(t,x,y,z);
      }
    }
  }
  return S;
}

const Field3D Wave1D::source_f(BoutReal t) {
  BoutReal x;
  Field3D result;
  result.allocate();

  int xi,yj,zk;

  for(xi=mesh->xstart;xi<mesh->xend+1;xi++)
    for(yj=mesh->ystart;yj < mesh->yend+1;yj++){
      for(zk=0;zk<mesh->LocalNz;zk++){
        x = mesh->GlobalX(xi)*Lx;
        result(xi,yj,zk) = -0.8*x*Cos(7*t)*Cos(2.0*Power(x,2)) - 2.0*Sin(10*t)*Sin(5.0*Power(x,2)) - 0.7;
      }
    }
  return result;
}

/////////////////// SOLUTION FOR G //////////////////////////////

//Manufactured solution
BoutReal MS_g(BoutReal t, BoutReal  x, BoutReal  UNUSED(y), BoutReal  UNUSED(z)) {
  // Input is in normalised x,y,z location
  x *= Lx;         // X input [0,1]
  return (BoutReal)0.7*x + 0.2*Sin(2.0*Power(x,2))*Cos(7*t) + 0.9;
}

//x-derivative of MS. For Neumann bnd cond
BoutReal dxMS_g(BoutReal t, BoutReal x, BoutReal  UNUSED(y), BoutReal  UNUSED(z)) {
  x *= Lx;         // X input [0,1]
  return 0.8*x*Cos(7*t)*Cos(2.0*Power(x,2)) + 0.7;
}

//Manufactured solution
const Field3D Wave1D::solution_g(BoutReal t) {
  Field3D S;
  S.allocate();
  
  S.setLocation(CELL_XLOW);
  
  int bx = (mesh->LocalNx - (mesh->xend - mesh->xstart + 1)) / 2;
  int by = (mesh->LocalNy - (mesh->yend - mesh->ystart + 1)) / 2;

  for (int xi = mesh->xstart - bx; xi < mesh->xend + bx + 1; xi++){
    for (int yj = mesh->ystart - by; yj < mesh->yend + by + 1; yj++){
      BoutReal x = mesh->GlobalX(xi);
      if(mesh->StaggerGrids) {
	x = 0.5*(mesh->GlobalX(xi-1) + mesh->GlobalX(xi));
      }
      BoutReal y = mesh->GlobalY(yj);//GlobalY not fixed yet
      for (int zk = 0; zk < mesh->LocalNz; zk++) {
        BoutReal z = coord->dz(xi,yj,zk)*zk;
        S(xi, yj, zk) = MS_g(t,x,y,z);
      }
    }
  }
  return S;
}

const Field3D Wave1D::source_g(BoutReal t) {
  BoutReal x;
  Field3D result;
  result.allocate();

  result.setLocation(CELL_XLOW);

  int xi,yj,zk;

  for(xi=mesh->xstart;xi<mesh->xend+1;xi++)
    for(yj=mesh->ystart;yj < mesh->yend+1;yj++){
      for(zk=0;zk<mesh->LocalNz;zk++){
	if(mesh->StaggerGrids) {
	  x = 0.5*(mesh->GlobalX(xi-1)+mesh->GlobalX(xi))*Lx;
	}else {
	  x = mesh->GlobalX(xi)*Lx;
	}
        result(xi,yj,zk) = -2.0*x*Cos(10*t)*Cos(5.0*Power(x,2)) - 1.4*Sin(7*t)*Sin(2.0*Power(x,2)) - 0.9;
      }
    }
  return result;
}
 


///////////////////////////////////////////////////////

BOUTMAIN(Wave1D); // Create a main() function
