
#include "bout/physicsmodel.hxx"

class FCIwave : public PhysicsModel {
private:
  Field3D n, nv; //< Evolving density, momentum

  Field3D Bxyz; ///< Total magnetic field
  
  /// Parallel divergence, using integration over projected cells
  Field3D Div_par_integrate(const Field3D &f) {
    Field3D f_B = f / Bxyz;
    
    f_B.splitYupYdown();
    mesh->getParallelTransform().integrateYUpDown(f_B);
    //mesh->getParallelTransform().calcYUpDown(f_B);
    
    f_B.applyParallelBoundary("parallel_dirichlet");

    Field3D result;
    result.allocate();
    
    Coordinates *coord = mesh->coordinates();
    
    for(auto i : result.region(RGN_NOBNDRY)) {
      result[i] = Bxyz[i] * (f_B.yup()[i.yp()] - f_B.ydown()[i.ym()])
        / (2.*coord->dy[i] * sqrt(coord->g_22[i]));

      if (!finite(result[i])) {
        output.write("[%d,%d,%d]: %e, %e -> %e\n",
                     i.x, i.y, i.z,
                     f_B.yup()[i.yp()],
                     f_B.ydown()[i.ym()],
                     result[i]);
        
      }
    }

    return result;
  }
  
protected:
  int init(bool restarting) override {

    // Get the magnetic field
    mesh->get(Bxyz, "B");
    
    Bxyz.applyBoundary("neumann");
    SAVE_ONCE(Bxyz);
    
    SOLVE_FOR2(n,nv);

    SAVE_REPEAT2(ddt(n), ddt(nv));
    
    return 0;
  }
  
  int rhs(BoutReal t) override {
    mesh->communicate(n,nv);

    n.applyParallelBoundary();
    nv.applyParallelBoundary();
    
    // Momentum
    ddt(nv) =
      - Div_par_integrate(SQ(nv)/floor(n, 1e-4))
      - Grad_par(n)
      + Grad2_par2(nv)
      ;

    // Density
    //ddt(n) = Div_par(nv);
    ddt(n) = -Div_par_integrate(nv);
    
    return 0;
  }
};

BOUTMAIN(FCIwave);
