
#include "bout/physicsmodel.hxx"

class FCIwave : public PhysicsModel {
private:
  Field3D logn, v; //< Evolving density, momentum
  Field3D n;
  
  Field3D Bxyz; ///< Total magnetic field

  bool expand_divergence;
  BoutReal background; ///< background density floor
  BoutReal log_background; // Log(background)
  
  /// Parallel divergence, using integration over projected cells
  Field3D Div_par_integrate(const Field3D &f) {
    Field3D f_B = f / Bxyz;
    
    f_B.splitYupYdown();
    mesh->getParallelTransform().integrateParallelSlices(f_B);

    // integrateParallelSlices replaces all yup/down points, so the boundary conditions
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
    
    Coordinates *coord = mesh->getCoordinates();
    
    for(auto i : result.getRegion(RGN_NOBNDRY)) {
      result[i] = Bxyz[i] * (f_B.yup()[i.yp()] - f_B.ydown()[i.ym()])
        / (2.*coord->dy[i] * sqrt(coord->g_22[i]));

      if (!finite(result[i])) {
        output.write("[%d,%d,%d]: %e, %e -> %e\n",
                     i.x(), i.y(), i.z(),
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

    Options::getRoot()->getSection("fciwave")->get("expand_divergence", expand_divergence, true);
    Options::getRoot()->getSection("fciwave")->get("background", background, 1e-6);
    log_background = log(background);
    
    SOLVE_FOR2(logn,v);

    SAVE_REPEAT(n);
    
    SAVE_REPEAT2(ddt(logn), ddt(v));
    
    return 0;
  }
  
  int rhs(BoutReal t) override {
    mesh->communicate(logn,v);

    // Boundary condition applied to log(n) to prevent negative densities
    logn.applyParallelBoundary();
    v.applyParallelBoundary();

    n = exp(logn);
    
    // Momentum
    ddt(v) =
      - v*Grad_par(v)
      - Grad_par(logn)
      + Grad2_par2(v)
      ;

    if (expand_divergence) {
      // Split the divergence of flux into two terms 
      ddt(logn) =
        -v * Grad_par(logn) - Div_par(v);
      
    } else {
      // Calculate the flux divergence using Div_par_integrate
      
      Field3D nv = n * v;
      nv.splitYupYdown();
      for (const auto &reg : mesh->getBoundariesPar()) {
        Field3D &nv_next = nv.ynext(reg->dir);
        nv_next.allocate();
        
        const Field3D &logn_next = logn.ynext(reg->dir);
        const Field3D &v_next = v.ynext(reg->dir);
        
        for (reg->first(); !reg->isDone(); reg->next()) {
          BoutReal n_b = exp(0.5*(logn_next(reg->x, reg->y+reg->dir, reg->z) +
                                  logn(reg->x, reg->y, reg->z)));
          BoutReal v_b = 0.5*(v_next(reg->x, reg->y+reg->dir, reg->z) +
                              v(reg->x, reg->y, reg->z));
          
          nv_next(reg->x, reg->y+reg->dir, reg->z) = 
            2.*n_b*v_b - nv(reg->x, reg->y, reg->z);
        }
      }
      
      // Logarithm of density
      ddt(logn) =
        - Div_par_integrate(nv) / floor(n, background);
      
      // Apply a soft floor to the density
      // Hard floors (setting ddt = 0) can slow convergence of solver
      for (auto i : logn.getRegion(RGN_NOBNDRY)) {
        if (ddt(logn)[i] < 0.0) {
          ddt(logn)[i] *= (1. - exp(log_background - logn[i]));
        }
      }
    }

    
    return 0;
  }
};

BOUTMAIN(FCIwave);
