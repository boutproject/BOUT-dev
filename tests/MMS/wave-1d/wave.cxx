#include <bout/physicsmodel.hxx>
#include <initialprofiles.hxx>
#include <derivs.hxx>
#include <math.h>
#include "mathematica.h"
#include <bout/constants.hxx>

BoutReal MS_f(BoutReal t, BoutReal  x, BoutReal  y, BoutReal  z);
BoutReal dxMS_f(BoutReal t, BoutReal  x, BoutReal  y, BoutReal  z);

BoutReal MS_g(BoutReal t, BoutReal  x, BoutReal  y, BoutReal  z);
BoutReal dxMS_g(BoutReal t, BoutReal  x, BoutReal  y, BoutReal  z);

BoutReal Lx, Ly, Lz; // Size of the domain


const Field3D HLL(const Field3D &f, const Field3D &u, BoutReal SL, BoutReal SR) {
  // A two-wave approximate solver, given left and right wave speeds such that
  // SL <= 0 <= SR
  // If SR = -SL  is the CFL limiting fastest wave speed, then this is
  // equivalent to the local Lax-Friedrichs or Rusanov method
  // First-order accurate
  Field3D result;
  result.allocate();
  
  Coordinates *coord = mesh->coordinates();
  
  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=mesh->ystart; j<=mesh->yend; j++)
      for(int k=0;k<mesh->LocalNz;k++) {
          
        BoutReal fp = (SR*f(i,j,k) - SL*f(i+1,j,k) + SL*SR*(u(i+1,j,k) - u(i,j,k)) ) / (SR - SL);
          BoutReal fm = (SR*f(i-1,j,k) - SL*f(i,j,k) + SL*SR*(u(i,j,k) - u(i-1,j,k)) ) / (SR - SL);
          
          result(i,j,k) = (fm - fp) / coord->dx(i,j);
        }
  return result;
}

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
  int init(bool restarting) {
    // Coordinate system
    coord = mesh->coordinates();
    
    // Get the options
    Options *meshoptions = Options::getRoot()->getSection("mesh");
    
    meshoptions->get("Lx",Lx,1.0);
    meshoptions->get("Ly",Ly,1.0);
    
    /*this assumes equidistant grid*/
    int nguard = mesh->xstart;

    Field2D dx=Lx/(mesh->GlobalNx - 2*nguard);
    // Set dx for the CELL_CENTRE location
    coord->dx = dx;
    // As dx is constant, we can just change the location, without interp_to
    dx.setLocation(CELL_XLOW);
    // Set dx for the CELL_XLOW location
    // If this would not be set, interp_to would be called.
    // As we have only 1 guard cell, that could fail.
    coord->dx.set(dx);
    coord->dy = Ly/(mesh->GlobalNy - 2*nguard);
    
    SAVE_ONCE2(Lx,Ly);
    
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
    
    //Dirichlet everywhere except inner x-boundary Neumann
    f.addBndryFunction(MS_f,BNDRY_ALL);
    f.addBndryFunction(dxMS_f,BNDRY_XIN);

    g.addBndryFunction(MS_g,BNDRY_ALL);
    g.addBndryFunction(dxMS_g,BNDRY_XIN);
    
    // Tell BOUT++ to solve f and g
    bout_solve(f, "f");
    bout_solve(g, "g");

    //Set initial condition to MS at t = 0.
    if(mesh->StaggerGrids) {
      for (int xi = mesh->xstart; xi < mesh->xend +1; xi++){
	for (int yj = mesh->ystart; yj < mesh->yend + 1; yj++){
	  for (int zk = 0; zk < mesh->LocalNz; zk++) {
	    f(xi, yj, zk) = MS_f(0.,mesh->GlobalX(xi),mesh->GlobalY(yj),coord->dz*zk);
	    g(xi, yj, zk) = MS_g(0.,0.5*(mesh->GlobalX(xi)+mesh->GlobalX(xi-1)),mesh->GlobalY(yj),coord->dz*zk);
	  }
	}
      }
    }else {
      for (int xi = mesh->xstart; xi < mesh->xend +1; xi++){
	for (int yj = mesh->ystart; yj < mesh->yend + 1; yj++){
	  for (int zk = 0; zk < mesh->LocalNz; zk++) {
	    f(xi, yj, zk) = MS_f(0.,mesh->GlobalX(xi),mesh->GlobalY(yj),coord->dz*zk);
	    g(xi, yj, zk) = MS_g(0.,mesh->GlobalX(xi),mesh->GlobalY(yj),coord->dz*zk);
	  }
	}
      }
    }
    E_f.allocate();
    E_g.allocate();
    SAVE_REPEAT2(E_f, E_g);
    
    return 0;
  }
  
  int rhs(BoutReal t) {
    mesh->communicate(f,g); // Communicate guard cells
    
    //update time-dependent boundary conditions
    f.applyBoundary(t);
    g.applyBoundary(t);
    
    // HLL. Related to Lax-Friedrichs
    //ddt(f) = HLL(-g, f, -1.0, 1.0);
    //ddt(g) = HLL(-f, g, -1.0, 1.0);

    // Central differencing
    ddt(f) = DDX(g, CELL_CENTRE);// + 20*SQ(coord->dx)*D2DX2(f);
    ddt(g) = DDX(f, CELL_XLOW);// + 20*SQ(coord->dx)*D2DX2(g);
    
    //add MMS source term
    ddt(f) += source_f(t);
    ddt(g) += source_g(t);
    
    return 0;
  }
  
  // This called every output timestep
  int outputMonitor(BoutReal simtime, int iter, int NOUT) {

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
          output_error.write("Error at %d,%d,%d = %e, %e",
                       xi, yj, zk, E_f(xi,yj,zk), E_g(xi,yj,zk));
        }
      }
    }
    return 0;
  }
};


/////////////////// SOLUTION FOR F //////////////////////////////

//Manufactured solution
BoutReal MS_f(BoutReal t, BoutReal  x, BoutReal  y, BoutReal  z) {
  // Input is in normalised x,y,z location
  x *= Lx;         // X input [0,1]
  y *= Ly / TWOPI; // Y input [0, 2pi]
  return (BoutReal)0.9 + 0.9*x + 0.2*Cos(10*t)*Sin(5.*Power(x,2));
}

//x-derivative of MS. For Neumann bnd cond
BoutReal dxMS_f(BoutReal t, BoutReal  x, BoutReal  y, BoutReal  z) {
  x *= Lx;         // X input [0,1]
  y *= Ly / TWOPI; // Y input [0, 2pi]
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
        BoutReal z = coord->dz*zk;
        S(xi, yj, zk) = MS_f(t,x,y,z);
      }
    }
  }
  return S;
}

const Field3D Wave1D::source_f(BoutReal t) {
  BoutReal x,y,z;
  Field3D result;
  result.allocate();

  int xi,yj,zk;

  for(xi=mesh->xstart;xi<mesh->xend+1;xi++)
    for(yj=mesh->ystart;yj < mesh->yend+1;yj++){
      for(zk=0;zk<mesh->LocalNz;zk++){
        x = mesh->GlobalX(xi)*Lx;
        y = mesh->GlobalY(yj)*Ly;
        z = zk*coord->dz;
        result(xi,yj,zk) = -0.8*x*Cos(7*t)*Cos(2.0*Power(x,2)) - 2.0*Sin(10*t)*Sin(5.0*Power(x,2)) - 0.7;
      }
    }
  return result;
}

/////////////////// SOLUTION FOR G //////////////////////////////

//Manufactured solution
BoutReal MS_g(BoutReal t, BoutReal  x, BoutReal  y, BoutReal  z) {
  // Input is in normalised x,y,z location
  x *= Lx;         // X input [0,1]
  y *= Ly / TWOPI; // Y input [0, 2pi]
  return (BoutReal)0.7*x + 0.2*Sin(2.0*Power(x,2))*Cos(7*t) + 0.9;
}

//x-derivative of MS. For Neumann bnd cond
BoutReal dxMS_g(BoutReal t, BoutReal  x, BoutReal  y, BoutReal  z) {
  x *= Lx;         // X input [0,1]
  y *= Ly / TWOPI; // Y input [0, 2pi]
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
        BoutReal z = coord->dz*zk;
        S(xi, yj, zk) = MS_g(t,x,y,z);
      }
    }
  }
  return S;
}

const Field3D Wave1D::source_g(BoutReal t) {
  BoutReal x,y,z;
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
        y = mesh->GlobalY(yj)*Ly;
        z = zk*coord->dz;
        result(xi,yj,zk) = -2.0*x*Cos(10*t)*Cos(5.0*Power(x,2)) - 1.4*Sin(7*t)*Sin(2.0*Power(x,2)) - 0.9;
      }
    }
  return result;
}
 


///////////////////////////////////////////////////////

BOUTMAIN(Wave1D); // Create a main() function

