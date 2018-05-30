
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

    // integrateYUpDown replaces all yup/down points, so the boundary conditions
    // now need to be applied. If Bxyz has neumann parallel boundary conditions
    // then the boundary condition is simpler since f = 0 gives f_B=0 boundary condition.

    /// Loop over the mesh boundary regions
    for (const auto &reg : mesh->getBoundariesPar()) {
      Field3D &f_B_next = f_B.ynext(reg->dir);
      const Field3D &f_next = f.ynext(reg->dir);
      const Field3D &B_next = Bxyz.ynext(reg->dir);
      
      for (reg->first(); !reg->isDone(); reg->next()) {
        f_B_next(reg->x, reg->y+reg->dir, reg->z) =
          f_next(reg->x, reg->y+reg->dir, reg->z) / B_next(reg->x, reg->y+reg->dir, reg->z);
      }
    }
    
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

    // Neumann boundaries simplifies parallel derivatives
    Bxyz.applyBoundary("neumann");
    Bxyz.applyParallelBoundary("parallel_neumann");
    SAVE_ONCE(Bxyz);
    
    SOLVE_FOR2(n,nv);

    SAVE_REPEAT2(ddt(n), ddt(nv));
    
    return 0;
  }
  
  int rhs(BoutReal t) override {
    mesh->communicate(n,nv);

    n.applyParallelBoundary();
    nv.applyParallelBoundary();

    // Calculate momentum flux
    Field3D momflux = SQ(nv)/floor(n, 1e-4);
    momflux.splitYupYdown();
    momflux.yup().allocate();
    momflux.ydown().allocate();
    momflux.applyParallelBoundary("parallel_dirichlet");
    
    // Momentum
    ddt(nv) =
      - Div_par_integrate(momflux)
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
