
#include "bout/physicsmodel.hxx"

class FCIwave : public PhysicsModel {
private:
  Field3D logn, nv; //< Evolving density, momentum
  Field3D n, v;
  
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
    
    SOLVE_FOR2(logn,nv);

    SAVE_REPEAT2(n, v)
    
    SAVE_REPEAT2(ddt(logn), ddt(nv));
    
    v.setBoundary("v");
    
    return 0;
  }
  
  int rhs(BoutReal t) override {
    mesh->communicate(logn,nv);

    // Boundary condition applied to log(n) to prevent negative densities
    logn.applyParallelBoundary();
    nv.applyParallelBoundary();

    n = exp(logn);
    v = nv / n;

    mesh->communicate(v);
    // Apply boundary condition to velocity v
    v.applyParallelBoundary("parallel_dirichlet");
    
    // Calculate momentum flux
    Field3D momflux = nv * v;
    momflux.splitYupYdown();
    momflux.yup().allocate();
    momflux.ydown().allocate();
    momflux.applyParallelBoundary("parallel_dirichlet");
    
    // Momentum
    ddt(nv) =
      - Div_par_integrate(momflux)
      - n * Grad_par(logn)
      + Grad2_par2(nv)
      ;

    // Density
    //ddt(n) = Div_par(nv);
    ddt(logn) =
      //-Div_par_integrate(nv) / n;
      -v * Grad_par(logn) - Div_par_integrate(v);
    
    
    return 0;
  }
};

BOUTMAIN(FCIwave);
