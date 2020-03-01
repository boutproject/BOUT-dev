#include <bout.hxx>
#include <boutmain.hxx>
#include <initialprofiles.hxx>
#include <derivs.hxx>
#include <cmath>
#include "mathematica.h"
#include <bout/constants.hxx>
#include <unused.hxx>

void solution(Field3D &f, BoutReal t, BoutReal D);
class ErrorMonitor: public Monitor{
public:
  int call(Solver *solver, BoutReal simtime, int iter, int NOUT) override;
};
BoutReal MS(BoutReal t, BoutReal  x, BoutReal  y, BoutReal  z);
BoutReal dxMS(BoutReal t, BoutReal  x, BoutReal  y, BoutReal  z);
Field3D MMS_Source(BoutReal t);


Field3D N;
Field3D E_N, S,source; //N error vector, solution,source

BoutReal mu_N; // Parallel collisional diffusion coefficient
BoutReal Lx, Ly, Lz;

Coordinates *coord;
ErrorMonitor error_monitor;
int physics_init(bool UNUSED(restarting)) {
  // Get the options
  Options *meshoptions = Options::getRoot()->getSection("mesh");

  coord = mesh->getCoordinates();
  
  meshoptions->get("Lx",Lx,1.0);
  meshoptions->get("Ly",Ly,1.0);

  /*this assumes equidistant grid*/
  int nguard = mesh->xstart;
  coord->dx = Lx/(mesh->GlobalNx - 2*nguard);
  coord->dy = Ly/(mesh->GlobalNy - 2*nguard);

  SAVE_ONCE2(Lx,Ly);

  Options *cytooptions = Options::getRoot()->getSection("cyto");
  cytooptions->get("dis", mu_N, 1);

  SAVE_ONCE(mu_N);

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

  //Dirichlet everywhere except inner x-boundary Neumann
  N.addBndryFunction(MS,BNDRY_ALL);
  N.addBndryFunction(dxMS,BNDRY_XIN);

  //Dirichlet boundary conditions everywhere
  //N.addBndryFunction(MS,BNDRY_ALL);

  // Tell BOUT++ to solve N
  SOLVE_FOR(N);

  //Set initial condition to MS at t = 0.
  for (int xi = mesh->xstart; xi < mesh->xend +1; xi++){
    for (int yj = mesh->ystart; yj < mesh->yend + 1; yj++){
      for (int zk = 0; zk < mesh->LocalNz; zk++) {
        output.write("Initial condition at {:d},{:d},{:d}\n", xi, yj, zk);
        N(xi, yj, zk) = MS(0.,mesh->GlobalX(xi)*Lx,mesh->GlobalY(yj)*Ly,coord->dz*zk);
      }
    }
  }

  E_N.allocate();
  SAVE_REPEAT(E_N);
  S.allocate();
  SAVE_REPEAT(S);
  source.allocate();
  SAVE_REPEAT(source);

  error_monitor.call(nullptr, 0,  0, 0);
  solver->addMonitor(&error_monitor);

  return 0;
}

int physics_run(BoutReal t) {
  mesh->communicate(N); // Communicate guard cells

  //update time-dependent boundary conditions
  output.write("APPLYING BOUNDARY\n");
  N.applyBoundary(t);

  ddt(N) = mu_N* D2DX2(N);


  //add MMS source term
  ddt(N) += MMS_Source(t);
  output.write("FINISHED RHS\n");
  return 0;
}


//Manufactured solution
BoutReal MS(BoutReal t, BoutReal  x, BoutReal  y, BoutReal  z) {
  output.write("-> MS at {:e}, {:e}, {:e}, {:e}\n", t, x, y, z);
  // Input is in normalised x,y,z location
  x *= Lx;         // X input [0,1]
  y *= Ly / TWOPI; // Y input [0, 2pi]
  return (BoutReal)0.9 + 0.9*x + 0.2*Cos(10*t)*Sin(5.*Power(x,2));
}

//x-derivative of MS. For Neumann bnd cond
BoutReal dxMS(BoutReal t, BoutReal  x, BoutReal  y, BoutReal  z) {
  output.write("-> dxMS at {:e}, {:e}, {:e}, {:e}\n", t, x, y, z);
  x *= Lx;         // X input [0,1]
  y *= Ly / TWOPI; // Y input [0, 2pi]
  return 0.9 + 2.*x*Cos(10*t)*Cos(5.*Power(x,2));
}


//Manufactured solution
void solution(Field3D &f, BoutReal t, BoutReal UNUSED(D)) {
  int bx = (mesh->LocalNx - (mesh->xend - mesh->xstart + 1)) / 2;
  int by = (mesh->LocalNy - (mesh->yend - mesh->ystart + 1)) / 2;
  BoutReal x,y,z;

  for (int xi = mesh->xstart - bx; xi < mesh->xend + bx + 1; xi++){
    for (int yj = mesh->ystart - by; yj < mesh->yend + by + 1; yj++){
      x = mesh->GlobalX(xi);
      y = mesh->GlobalY(yj);//GlobalY not fixed yet
      for (int zk = 0; zk < mesh->LocalNz; zk++) {
        z = coord->dz*zk;
        output.write("Solution at {:d},{:d},{:d}\n", xi, yj, zk);
        f(xi, yj, zk) = MS(t,x,y,z);
      }
    }

  }

}

//Source term calculated in and cut and pasted from Mathematica.
//\partial_t MS - \mu_N \partial^2_{xx } N = MMS_Source
Field3D MMS_Source(BoutReal t)
{
  BoutReal x;
  Field3D result;
  result = 0.0;
  
  int xi,yj,zk;

  for(xi=mesh->xstart;xi<mesh->xend+1;xi++)
    for(yj=mesh->ystart;yj < mesh->yend+1;yj++){
      for(zk=0;zk<mesh->LocalNz;zk++){
        x = mesh->GlobalX(xi)*Lx;
        result(xi,yj,zk) = -2.*Sin(10*t)*Sin(5.*Power(x,2)) + Cos(10*t)*
          (-2.*Cos(5.*Power(x,2)) + 20.*Power(x,2)*Sin(5.*Power(x,2)));
      }
    }
  return result;
}

int ErrorMonitor::call(Solver *UNUSED(solver), BoutReal simtime, int UNUSED(iter), int UNUSED(NOUT)) {
  solution(S, simtime, mu_N);

  //Calculate the error. norms are calculated in the post-processing
  E_N = 0.;
  for (int xi = mesh->xstart; xi < mesh->xend + 1; xi++){
    for (int yj = mesh->ystart ; yj < mesh->yend + 1; yj++){
      for (int zk = 0; zk < mesh->LocalNz ; zk++) {
        E_N(xi, yj, zk) = N(xi, yj, zk) - S(xi, yj, zk);

        output_error.write("Error({:d},{:d},{:d}): {:e}, {:e} -> {:e}\n",
                     xi, yj, zk, 
                     N(xi, yj, zk), S(xi, yj, zk), E_N(xi, yj, zk));
      }
    }
  }
  source = MMS_Source(simtime);
  return 0;
}

