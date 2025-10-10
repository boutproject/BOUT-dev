
#include "bout/parallel_boundary_region.hxx"
#include "bout/physicsmodel.hxx"

class FCIwave : public PhysicsModel {
private:
  Field3D n, nv; //< Evolving density, momentum
  Field3D logn, v;

  Field3D Bxyz; ///< Total magnetic field

  bool div_integrate;      ///< Use area integration for divergence operator in density
  bool log_density;        ///< Evolve logarithm of density
  BoutReal background;     ///< background density floor
  BoutReal log_background; // Log(background)

  /// Parallel divergence, using integration over projected cells
  Field3D Div_par_integrate(const Field3D& f) {
    Field3D f_B = f / Bxyz;

    f_B.splitParallelSlices();
    mesh->getCoordinates()->getParallelTransform().integrateParallelSlices(f_B);

    // integrateParallelSlices replaces all yup/down points, so the boundary conditions
    // now need to be applied. If Bxyz has neumann parallel boundary conditions
    // then the boundary condition is simpler since f = 0 gives f_B=0 boundary condition.

    /// Loop over the mesh boundary regions
    for (const auto& reg : mesh->getBoundariesPar()) {
      Field3D& f_B_next = f_B.ynext(reg->dir);
      const Field3D& f_next = f.ynext(reg->dir);
      const Field3D& B_next = Bxyz.ynext(reg->dir);

      for (reg->first(); !reg->isDone(); reg->next()) {
        f_B_next(reg->ind().x(), reg->ind().y() + reg->dir, reg->ind().z()) =
            f_next(reg->ind().x(), reg->ind().y() + reg->dir, reg->ind().z())
            / B_next(reg->ind().x(), reg->ind().y() + reg->dir, reg->ind().z());
      }
    }

    Field3D result;
    result.allocate();

    Coordinates* coord = mesh->getCoordinates();

    for (auto i : result.getRegion(RGN_NOBNDRY)) {
      result[i] = Bxyz[i] * (f_B.yup()[i.yp()] - f_B.ydown()[i.ym()])
                  / (2. * coord->dy[i] * sqrt(coord->g_22[i]));

      if (!finite(result[i])) {
        output.write("[{:d},{:d},{:d}]: {:e}, {:e} -> {:e}\n", i.x(), i.y(), i.z(),
                     f_B.yup()[i.yp()], f_B.ydown()[i.ym()], result[i]);
      }
    }

    return result;
  }

protected:
  int init(bool UNUSED(restarting)) override {

    // Get the magnetic field
    mesh->get(Bxyz, "B");

    auto& options = Options::root()["fciwave"];
    div_integrate = options["div_integrate"].withDefault(true);
    log_density = options["log_density"].withDefault(false);
    background = options["background"].withDefault(false);
    log_background = log(background);

    // Neumann boundaries simplifies parallel derivatives
    Bxyz.applyBoundary("neumann");
    Bxyz.applyParallelBoundary("parallel_neumann_o2");
    SAVE_ONCE(Bxyz);

    SOLVE_FOR(nv);
    if (log_density) {
      SOLVE_FOR(logn);
      SAVE_REPEAT(n);
    } else {
      SOLVE_FOR(n);
    }

    v.setBoundary("v");

    return 0;
  }

  int rhs(BoutReal UNUSED(time)) override {
    if (log_density) {
      mesh->communicate(logn, nv);
      // Apply boundary condition to log(n)
      // rather than n to prevent negative densities
      logn.applyParallelBoundary();

      n = exp(logn);
      n.splitParallelSlices();
      n.yup() = exp(logn.yup());
      n.ydown() = exp(logn.ydown());
    } else {
      mesh->communicate(n, nv);

      n.applyParallelBoundary();
    }

    // Calculate velocity and momentum flux
    v = nv / floor(n, 1e-4);
    Field3D momflux = nv * v;

    // Apply boundary conditions to v
    v.splitParallelSlices();
    v.yup().allocate();
    v.ydown().allocate();
    v.applyParallelBoundary();

    // Ensure that boundary conditions are consistent
    // between v, nv and momentum flux

    momflux.splitParallelSlices();
    for (const auto& reg : mesh->getBoundariesPar()) {
      // Using the values of density and velocity on the boundary
      const Field3D& n_next = n.ynext(reg->dir);
      const Field3D& v_next = v.ynext(reg->dir);

      // Set the momentum and momentum flux
      Field3D& nv_next = nv.ynext(reg->dir);
      Field3D& momflux_next = momflux.ynext(reg->dir);
      momflux_next.allocate();

      for (reg->first(); !reg->isDone(); reg->next()) {
        // Density at the boundary
        // Note: If evolving density, this should interpolate logn
        // but neumann boundaries are used here anyway.
        BoutReal n_b =
            0.5
            * (n_next(reg->ind().x(), reg->ind().y() + reg->dir, reg->ind().z())
               + n(reg->ind().x(), reg->ind().y(), reg->ind().z()));
        // Velocity at the boundary
        BoutReal v_b =
            0.5
            * (v_next(reg->ind().x(), reg->ind().y() + reg->dir, reg->ind().z())
               + v(reg->ind().x(), reg->ind().y(), reg->ind().z()));

        nv_next(reg->ind().x(), reg->ind().y() + reg->dir, reg->ind().z()) =
            2. * n_b * v_b - nv(reg->ind().x(), reg->ind().y(), reg->ind().z());

        momflux_next(reg->ind().x(), reg->ind().y() + reg->dir, reg->ind().z()) =
            2. * n_b * v_b * v_b
            - momflux(reg->ind().x(), reg->ind().y(), reg->ind().z());
      }
    }

    // Momentum
    ddt(nv) = -Div_par_integrate(momflux) - Grad_par(n) + Grad2_par2(nv);

    // Density
    if (div_integrate) {
      ddt(n) = -Div_par_integrate(nv);
    } else {
      ddt(n) = -Div_par(nv);
    }

    if (log_density) {
      ddt(logn) = ddt(n) / floor(n, background);

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
