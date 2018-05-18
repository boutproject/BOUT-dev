
#include "bout/physicsmodel.hxx"



class FCItest : public PhysicsModel {
private:
  Field3D f, g;

  Field3D Bxyz; ///< Total magnetic field

  /// Parallel divergence, using integration over projected cells
  Field3D Div_par_integrate(const Field3D &f) {
    Field3D f_B = f / Bxyz;

    mesh->getParallelTransform().integrateYUpDown(f_B);

    Field3D result;
    result.allocate();

    Coordinates *coord = mesh->coordinates();
    
    for(auto i : result.region(RGN_NOBNDRY)) {
      result[i] = Bxyz[i] * (f_B.yup()[i.yp()] - f_B.ydown()[i.ym()])
        / (coord->dy[i] * sqrt(coord->g_22[i]));
    }

    return result;
  }
protected:
  int init(bool restarting) override {

    // Get the magnetic field
    mesh->get(Bxyz, "B");
    Bxyz.applyBoundary("neumann");
    SAVE_ONCE(Bxyz);
    
    SOLVE_FOR2(f,g);
    return 0;
  }
  
  int rhs(BoutReal t) override {
    mesh->communicate(f,g);
    
    ddt(f) = Grad_par(g);
    ddt(g) = Div_par_integrate(f);
    
    return 0;
  }
};

BOUTMAIN(FCItest);
