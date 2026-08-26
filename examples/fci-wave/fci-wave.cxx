#include <bout/boundary_region_iter.hxx>
#include <bout/bout_types.hxx>
#include <bout/coordinates.hxx>
#include <bout/difops.hxx>
#include <bout/field3d.hxx>
#include <bout/options.hxx>
#include <bout/output.hxx>
#include <bout/physicsmodel.hxx>
#include <bout/unused.hxx>
#include <bout/yboundary_regions.hxx>

#include <cmath>

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
    const Coordinates* coord = mesh->getCoordinates();
    coord->getParallelTransform().integrateParallelSlices(f_B);

    // integrateParallelSlices replaces all yup/down points, so the boundary conditions
    // now need to be applied. If Bxyz has neumann parallel boundary conditions
    // then the boundary condition is simpler since f = 0 gives f_B=0 boundary condition.

    const auto ybndry = coord->getYBoundary();
    ybndry.iter([&](const bout::boundary::BoundaryIterator auto& point) {
      point.next(f_B) = point.next(f) / point.next(Bxyz);
    });

    Field3D result;
    result.allocate();

    for (auto i : result.getRegion(RGN_NOBNDRY)) {
      result[i] = Bxyz[i] * (f_B.yup()[i.yp()] - f_B.ydown()[i.ym()])
                  / (2. * coord->dy()[i] * sqrt(coord->g_22()[i]));

      if (!std::isfinite(result[i])) {
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
    background = options["background"].withDefault(0.0);
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
    const auto ybndry = mesh->getCoordinates()->getYBoundary();
    ybndry.iter([&](const bout::boundary::BoundaryIterator auto& point) {
      // Using the values of density and velocity on the boundary
      const BoutReal n_b = 0.5 * point.next(n) + point.current(n);
      const BoutReal v_b = 0.5 * point.next(v) + point.current(v);

      // Set the momentum and momentum flux
      point.next(nv) = 2. * n_b * v_b - point.current(nv);
      point.next(momflux) = 2. * n_b * v_b * v_b - point.current(momflux);
    });

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
