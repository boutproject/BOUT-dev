#include <bout/fv_ops.hxx>
#include <bout/globals.hxx>
#include <bout/msg_stack.hxx>
#include <bout/output.hxx>
#include <bout/utils.hxx>

namespace {
template <class T>
struct Slices {
  T c;
  T up;
  T down;
  Slices(bool use_slices, const T& field)
      : c(use_slices ? field : toFieldAligned(field)), up(use_slices ? field.yup() : c),
        down(use_slices ? field.ydown() : c) {}
};
template <class T>
Slices<T> makeslices(bool use_slices, const T& field) {
  return Slices<T>(use_slices, field);
}
} // namespace

namespace FV {

// Div ( a Grad_perp(f) ) -- ∇ ⋅ ( a ∇⊥ f) --  Vorticity
Field3D Div_a_Grad_perp(const Field3D& a, const Field3D& f) {
  ASSERT2(a.getLocation() == f.getLocation());

  Mesh* mesh = a.getMesh();

  Field3D result{zeroFrom(f)};

  Coordinates* coord = f.getCoordinates();

  // Flux in x

  int xs = mesh->xstart - 1;
  int xe = mesh->xend;

  /*
    if(mesh->firstX())
    xs += 1;
  */
  /*
    if(mesh->lastX())
    xe -= 1;
  */

  for (int i = xs; i <= xe; i++) {
    for (int j = mesh->ystart; j <= mesh->yend; j++) {
      for (int k = 0; k < mesh->LocalNz; k++) {
        // Calculate flux from i to i+1

        BoutReal fout = 0.5 * (a(i, j, k) + a(i + 1, j, k))
                        * (coord->J(i, j, k) * coord->g11(i, j, k)
                           + coord->J(i + 1, j, k) * coord->g11(i + 1, j, k))
                        * (f(i + 1, j, k) - f(i, j, k))
                        / (coord->dx(i, j, k) + coord->dx(i + 1, j, k));

        result(i, j, k) += fout / (coord->dx(i, j, k) * coord->J(i, j, k));
        result(i + 1, j, k) -= fout / (coord->dx(i + 1, j, k) * coord->J(i + 1, j, k));
      }
    }
  }

  const bool fci = f.hasParallelSlices() && a.hasParallelSlices();

  if (bout::build::use_metric_3d and fci) {
    // 3D Metric, need yup/ydown fields.
    // Requires previous communication of metrics
    // -- should insert communication here?
    if (!coord->g23.hasParallelSlices() || !coord->g_23.hasParallelSlices()
        || !coord->dy.hasParallelSlices() || !coord->dz.hasParallelSlices()
        || !coord->Bxy.hasParallelSlices() || !coord->J.hasParallelSlices()) {
      throw BoutException("metrics have no yup/down: Maybe communicate in init?");
    }
  }

  // Y and Z fluxes require Y derivatives

  // Fields containing values along the magnetic field

  // Values on this y slice (centre).
  // This is needed because toFieldAligned may modify the field
  const auto f_slice = makeslices(fci, f);
  const auto a_slice = makeslices(fci, a);

  // Only in 3D case with FCI do the metrics have parallel slices
  const bool metric_fci = fci and bout::build::use_metric_3d;
  const auto g23 = makeslices(metric_fci, coord->g23);
  const auto g_23 = makeslices(metric_fci, coord->g_23);
  const auto J = makeslices(metric_fci, coord->J);
  const auto dy = makeslices(metric_fci, coord->dy);
  const auto dz = makeslices(metric_fci, coord->dz);
  const auto Bxy = makeslices(metric_fci, coord->Bxy);

  // Result of the Y and Z fluxes
  Field3D yzresult(0.0, mesh);
  if (!fci) {
    yzresult.setDirectionY(YDirectionType::Aligned);
  }

  // Y flux

  BOUT_FOR(i, yzresult.getRegion("RGN_NOBNDRY")) {
    auto ikp = i.zp();
    auto ikm = i.zm();
    {
      const BoutReal coef = 0.5
                            * (g_23.c[i] / SQ(J.c[i] * Bxy.c[i])
                               + g_23.up[i.yp()] / SQ(J.up[i.yp()] * Bxy.up[i.yp()]));

      // Calculate Z derivative at y boundary
      const BoutReal dfdz = 0.5
                            * (f_slice.c[ikp] - f_slice.c[ikm] + f_slice.up[ikp.yp()]
                               - f_slice.up[ikm.yp()])
                            / (dz.c[i] + dz.up[i.yp()]);

      // Y derivative
      const BoutReal dfdy =
          2. * (f_slice.up[i.yp()] - f_slice.c[i]) / (dy.up[i.yp()] + dy.c[i]);

      const BoutReal fout = 0.25 * (a_slice.c[i] + a_slice.up[i.yp()])
                            * (J.c[i] * g23.c[i] + J.up[i.yp()] * g23.up[i.yp()])
                            * (dfdz - coef * dfdy);

      yzresult[i] += fout / (dy.c[i] * J.c[i]);
    }
    {
      // Calculate flux between j and j-1
      const BoutReal coef =
          0.5
          * (g_23.c[i] / SQ(J.c[i] * Bxy.c[i])
             + g_23.down[i.ym()] / SQ(J.down[i.ym()] * Bxy.down[i.ym()]));

      const BoutReal dfdz = 0.5
                            * (f_slice.c[ikp] - f_slice.c[ikm] + f_slice.down[ikp.ym()]
                               - f_slice.down[ikm.ym()])
                            / (dz.c[i] + dz.down[i.ym()]);

      const BoutReal dfdy =
          2. * (f_slice.c[i] - f_slice.down[i.ym()]) / (dy.c[i] + dy.down[i.ym()]);

      const BoutReal fout = 0.25 * (a_slice.c[i] + a_slice.down[i.ym()])
                            * (J.c[i] * g23.c[i] + J.down[i.ym()] * g23.down[i.ym()])
                            * (dfdz - coef * dfdy);

      yzresult[i] -= fout / (dy.c[i] * J.c[i]);
    }
    // Z flux
    {

      // Coefficient in front of df/dy term
      const BoutReal coef = g_23.c[i] / (dy.up[i.yp()] + 2. * dy.c[i] + dy.down[i.ym()])
                            / SQ(J.c[i] * Bxy.c[i]);

      const BoutReal fout =
          0.25 * (a_slice.c[i] + a_slice.c[ikp])
          * (J.c[i] * coord->g33[i] + J.c[ikp] * coord->g33[ikp])
          * ( // df/dz
              (f_slice.c[ikp] - f_slice.c[i]) / dz.c[i]
              // - g_yz * df/dy / SQ(J*B)
              - coef
                    * (f_slice.up[i.yp()] + f_slice.up[ikp.yp()] - f_slice.down[i.ym()]
                       - f_slice.down[ikp.ym()]));

      yzresult[i] += fout / (J.c[i] * dz.c[i]);
      yzresult[ikp] -= fout / (J.c[ikp] * dz.c[ikp]);
    }
  }

  // Check if we need to transform back
  if (fci) {
    result += yzresult;
  } else {
    result += fromFieldAligned(yzresult);
  }

  return result;
}

const Field3D Div_par_K_Grad_par(const Field3D& Kin, const Field3D& fin,
                                 bool bndry_flux) {
  TRACE("FV::Div_par_K_Grad_par");

  ASSERT2(Kin.getLocation() == fin.getLocation());

  Mesh* mesh = Kin.getMesh();

  bool use_parallel_slices = (Kin.hasParallelSlices() && fin.hasParallelSlices());

  const auto& K = use_parallel_slices ? Kin : toFieldAligned(Kin, "RGN_NOX");
  const auto& f = use_parallel_slices ? fin : toFieldAligned(fin, "RGN_NOX");

  Field3D result{zeroFrom(f)};

  // K and f fields in yup and ydown directions
  const auto& Kup = use_parallel_slices ? Kin.yup() : K;
  const auto& Kdown = use_parallel_slices ? Kin.ydown() : K;
  const auto& fup = use_parallel_slices ? fin.yup() : f;
  const auto& fdown = use_parallel_slices ? fin.ydown() : f;

  Coordinates* coord = fin.getCoordinates();

  BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
    // Calculate flux at upper surface

    const auto iyp = i.yp();
    const auto iym = i.ym();

    if (bndry_flux || mesh->periodicY(i.x()) || !mesh->lastY(i.x())
        || (i.y() != mesh->yend)) {

      BoutReal c = 0.5 * (K[i] + Kup[iyp]);             // K at the upper boundary
      BoutReal J = 0.5 * (coord->J[i] + coord->J[iyp]); // Jacobian at boundary
      BoutReal g_22 = 0.5 * (coord->g_22[i] + coord->g_22[iyp]);

      BoutReal gradient = 2. * (fup[iyp] - f[i]) / (coord->dy[i] + coord->dy[iyp]);

      BoutReal flux = c * J * gradient / g_22;

      result[i] += flux / (coord->dy[i] * coord->J[i]);
    }

    // Calculate flux at lower surface
    if (bndry_flux || mesh->periodicY(i.x()) || !mesh->firstY(i.x())
        || (i.y() != mesh->ystart)) {
      BoutReal c = 0.5 * (K[i] + Kdown[iym]);           // K at the lower boundary
      BoutReal J = 0.5 * (coord->J[i] + coord->J[iym]); // Jacobian at boundary

      BoutReal g_22 = 0.5 * (coord->g_22[i] + coord->g_22[iym]);

      BoutReal gradient = 2. * (f[i] - fdown[iym]) / (coord->dy[i] + coord->dy[iym]);

      BoutReal flux = c * J * gradient / g_22;

      result[i] -= flux / (coord->dy[i] * coord->J[i]);
    }
  }

  if (!use_parallel_slices) {
    // Shifted to field aligned coordinates, so need to shift back
    result = fromFieldAligned(result, "RGN_NOBNDRY");
  }

  return result;
}

const Field3D D4DY4(const Field3D& d_in, const Field3D& f_in) {
  ASSERT1_FIELDS_COMPATIBLE(d_in, f_in);

  Mesh* mesh = d_in.getMesh();

  Coordinates* coord = f_in.getCoordinates();

  ASSERT2(d_in.getDirectionY() == f_in.getDirectionY());
  const bool are_unaligned = ((d_in.getDirectionY() == YDirectionType::Standard)
                              and (f_in.getDirectionY() == YDirectionType::Standard));

  // Convert to field aligned coordinates
  Field3D d = are_unaligned ? toFieldAligned(d_in, "RGN_NOX") : d_in;
  Field3D f = are_unaligned ? toFieldAligned(f_in, "RGN_NOX") : f_in;

  Field3D result{zeroFrom(f)};

  for (int i = mesh->xstart; i <= mesh->xend; i++) {
    // Check for boundaries
    bool yperiodic = mesh->periodicY(i);
    bool has_upper_boundary = !yperiodic && mesh->lastY(i);
    bool has_lower_boundary = !yperiodic && mesh->firstY(i);

    // Always calculate fluxes at upper Y cell boundary
    const int ystart =
        has_lower_boundary
            ? mesh->ystart
            : // Don't calculate flux from boundary mesh->ystart-1 into domain
            mesh->ystart - 1; // Calculate flux from last guard cell into domain

    const int yend = has_upper_boundary
                         ? mesh->yend - 1
                         : // Don't calculate flux from mesh->yend into boundary
                         mesh->yend;

    for (int j = ystart; j <= yend; j++) {
      for (int k = 0; k < mesh->LocalNz; k++) {
        BoutReal dy3 = SQ(coord->dy(i, j, k)) * coord->dy(i, j, k);
        // 3rd derivative at upper boundary

        BoutReal d3fdy3 =
            (f(i, j + 2, k) - 3. * f(i, j + 1, k) + 3. * f(i, j, k) - f(i, j - 1, k))
            / dy3;

        BoutReal flux = 0.5 * (d(i, j, k) + d(i, j + 1, k))
                        * (coord->J(i, j, k) + coord->J(i, j + 1, k)) * d3fdy3;

        result(i, j, k) += flux / (coord->J(i, j, k) * coord->dy(i, j, k));
        result(i, j + 1, k) -= flux / (coord->J(i, j + 1, k) * coord->dy(i, j + 1, k));
      }
    }
  }

  // Convert result back to non-aligned coordinates
  return are_unaligned ? fromFieldAligned(result, "RGN_NOBNDRY") : result;
}

const Field3D D4DY4_Index(const Field3D& f_in, bool bndry_flux) {
  Mesh* mesh = f_in.getMesh();

  // Convert to field aligned coordinates
  const bool is_unaligned = (f_in.getDirectionY() == YDirectionType::Standard);
  Field3D f = is_unaligned ? toFieldAligned(f_in, "RGN_NOX") : f_in;

  Field3D result{zeroFrom(f)};

  Coordinates* coord = f_in.getCoordinates();

  for (int i = mesh->xstart; i <= mesh->xend; i++) {
    bool yperiodic = mesh->periodicY(i);

    bool has_upper_boundary = !yperiodic && mesh->lastY(i);
    bool has_lower_boundary = !yperiodic && mesh->firstY(i);

    for (int j = mesh->ystart; j <= mesh->yend; j++) {

      // Right boundary

      if (bndry_flux || (j != mesh->yend) || !has_upper_boundary) {
        // Calculate the fluxes

        if (j != mesh->yend || !has_upper_boundary) {

          for (int k = 0; k < mesh->LocalNz; k++) {
            // Right boundary common factors
            const BoutReal common_factor = 0.25
                                           * (coord->dy(i, j, k) + coord->dy(i, j + 1, k))
                                           * (coord->J(i, j, j) + coord->J(i, j + 1, k));

            const BoutReal factor_rc =
                common_factor / (coord->J(i, j, k) * coord->dy(i, j, k));
            const BoutReal factor_rp =
                common_factor / (coord->J(i, j + 1, k) * coord->dy(i, j + 1, k));

            // Not on domain boundary
            // 3rd derivative at right cell boundary

            const BoutReal d3fdx3 =
                (f(i, j + 2, k) - 3. * f(i, j + 1, k) + 3. * f(i, j, k) - f(i, j - 1, k));

            result(i, j, k) += d3fdx3 * factor_rc;
            result(i, j + 1, k) -= d3fdx3 * factor_rp;
          }
        } else {
          // At a domain boundary
          // Use a one-sided difference formula

          for (int k = 0; k < mesh->LocalNz; k++) {
            // Right boundary common factors
            const BoutReal common_factor = 0.25
                                           * (coord->dy(i, j, k) + coord->dy(i, j + 1, k))
                                           * (coord->J(i, j, j) + coord->J(i, j + 1, k));

            const BoutReal factor_rc =
                common_factor / (coord->J(i, j, k) * coord->dy(i, j, k));
            const BoutReal factor_rp =
                common_factor / (coord->J(i, j + 1, k) * coord->dy(i, j + 1, k));

            const BoutReal d3fdx3 =
                -((16. / 5) * 0.5 * (f(i, j + 1, k) + f(i, j, k)) // Boundary value f_b
                  - 6. * f(i, j, k)                               // f_0
                  + 4. * f(i, j - 1, k)                           // f_1
                  - (6. / 5) * f(i, j - 2, k)                     // f_2
                );

            result(i, j, k) += d3fdx3 * factor_rc;
            result(i, j + 1, k) -= d3fdx3 * factor_rp;
          }
        }
      }

      // Left cell boundary

      if (bndry_flux || (j != mesh->ystart) || !has_lower_boundary) {
        // Calculate the fluxes

        if (j != mesh->ystart || !has_lower_boundary) {
          for (int k = 0; k < mesh->LocalNz; k++) {
            const BoutReal common_factor = 0.25
                                           * (coord->dy(i, j, k) + coord->dy(i, j + 1, k))
                                           * (coord->J(i, j, k) + coord->J(i, j - 1, k));

            const BoutReal factor_lc =
                common_factor / (coord->J(i, j, k) * coord->dy(i, j, k));
            const BoutReal factor_lm =
                common_factor / (coord->J(i, j - 1, k) * coord->dy(i, j - 1, k));

            // Not on a domain boundary
            const BoutReal d3fdx3 =
                (f(i, j + 1, k) - 3. * f(i, j, k) + 3. * f(i, j - 1, k) - f(i, j - 2, k));

            result(i, j, k) -= d3fdx3 * factor_lc;
            result(i, j - 1, k) += d3fdx3 * factor_lm;
          }
        } else {
          // On a domain (Y) boundary
          for (int k = 0; k < mesh->LocalNz; k++) {
            const BoutReal common_factor = 0.25
                                           * (coord->dy(i, j, k) + coord->dy(i, j + 1, k))
                                           * (coord->J(i, j, k) + coord->J(i, j - 1, k));

            const BoutReal factor_lc =
                common_factor / (coord->J(i, j, k) * coord->dy(i, j, k));
            const BoutReal factor_lm =
                common_factor / (coord->J(i, j - 1, k) * coord->dy(i, j - 1, k));
            const BoutReal d3fdx3 =
                -(-(16. / 5) * 0.5 * (f(i, j - 1, k) + f(i, j, k)) // Boundary value f_b
                  + 6. * f(i, j, k)                                // f_0
                  - 4. * f(i, j + 1, k)                            // f_1
                  + (6. / 5) * f(i, j + 2, k)                      // f_2
                );

            result(i, j, k) -= d3fdx3 * factor_lc;
            result(i, j - 1, k) += d3fdx3 * factor_lm;
          }
        }
      }
    }
  }

  // Convert result back to non-aligned coordinates
  return is_unaligned ? fromFieldAligned(result, "RGN_NOBNDRY") : result;
}

void communicateFluxes(Field3D& f) {
  Mesh* mesh = f.getMesh();

  // Use X=0 as temporary buffer
  if (mesh->xstart != 2) {
    throw BoutException("communicateFluxes: Sorry!");
  }

  int size = mesh->LocalNy * mesh->LocalNz;
  comm_handle xin, xout;
  // Cache results to silence spurious compiler warning about xin,
  // xout possibly being uninitialised when used
  const bool not_first = mesh->periodicX || !mesh->firstX();
  const bool not_last = mesh->periodicX || !mesh->lastX();

  if (not_first) {
    xin = mesh->irecvXIn(f(0, 0), size, 0);
  }
  if (not_last) {
    xout = mesh->irecvXOut(f(mesh->LocalNx - 1, 0), size, 1);
  }
  // Send X=1 values
  if (not_first) {
    mesh->sendXIn(f(1, 0), size, 1);
  }
  if (not_last) {
    mesh->sendXOut(f(mesh->LocalNx - 2, 0), size, 0);
  }
  // Wait
  if (not_first) {
    mesh->wait(xin);
    // Add to cells
    for (int y = mesh->ystart; y <= mesh->yend; y++) {
      for (int z = 0; z < mesh->LocalNz; z++) {
        f(2, y, z) += f(0, y, z);
      }
    }
  }
  if (not_last) {
    mesh->wait(xout);
    // Add to cells
    for (int y = mesh->ystart; y <= mesh->yend; y++) {
      for (int z = 0; z < mesh->LocalNz; z++) {
        f(mesh->LocalNx - 3, y, z) += f(mesh->LocalNx - 1, y, z);
      }
    }
  }
}

Field3D Div_Perp_Lap(const Field3D& a, const Field3D& f, CELL_LOC outloc) {

  Field3D result = 0.0;

  //////////////////////////////////////////
  // X-Z diffusion
  //
  //            Z
  //            |
  //
  //     o --- gU --- o
  //     |     nU     |
  //     |            |
  //    gL nL      nR gR    -> X
  //     |            |
  //     |     nD     |
  //     o --- gD --- o
  //
  Coordinates* coords = a.getCoordinates(outloc);
  Mesh* mesh = f.getMesh();

  for (int i = mesh->xstart; i <= mesh->xend; i++) {
    for (int j = mesh->ystart; j <= mesh->yend; j++) {
      for (int k = 0; k < mesh->LocalNz; k++) {

        // wrap k-index around as Z is (currently) periodic.
        int kp = (k + 1) % (mesh->LocalNz);
        int km = (k - 1 + mesh->LocalNz) % (mesh->LocalNz);

        // Calculate gradients on cell faces -- assumes constant grid spacing

        BoutReal gR =
            (coords->g11(i, j, k) + coords->g11(i + 1, j, k))
                * (f(i + 1, j, k) - f(i, j, k))
                / (coords->dx(i + 1, j, k) + coords->dx(i, j, k))
            + 0.5 * (coords->g13(i, j, k) + coords->g13(i + 1, j, k))
                  * (f(i + 1, j, kp) - f(i + 1, j, km) + f(i, j, kp) - f(i, j, km))
                  / (4. * coords->dz(i, j, k));

        BoutReal gL =
            (coords->g11(i - 1, j, k) + coords->g11(i, j, k))
                * (f(i, j, k) - f(i - 1, j, k))
                / (coords->dx(i - 1, j, k) + coords->dx(i, j, k))
            + 0.5 * (coords->g13(i - 1, j, k) + coords->g13(i, j, k))
                  * (f(i - 1, j, kp) - f(i - 1, j, km) + f(i, j, kp) - f(i, j, km))
                  / (4 * coords->dz(i, j, k));

        BoutReal gD =
            coords->g13(i, j, k)
                * (f(i + 1, j, km) - f(i - 1, j, km) + f(i + 1, j, k) - f(i - 1, j, k))
                / (4. * coords->dx(i, j, k))
            + coords->g33(i, j, k) * (f(i, j, k) - f(i, j, km)) / coords->dz(i, j, k);

        BoutReal gU =
            coords->g13(i, j, k)
                * (f(i + 1, j, kp) - f(i - 1, j, kp) + f(i + 1, j, k) - f(i - 1, j, k))
                / (4. * coords->dx(i, j, k))
            + coords->g33(i, j, k) * (f(i, j, kp) - f(i, j, k)) / coords->dz(i, j, k);

        // Flow right
        BoutReal flux = gR * 0.25 * (coords->J(i + 1, j, k) + coords->J(i, j, k))
                        * (a(i + 1, j, k) + a(i, j, k));
        result(i, j, k) += flux / (coords->dx(i, j, k) * coords->J(i, j, k));

        // Flow left
        flux = gL * 0.25 * (coords->J(i - 1, j, k) + coords->J(i, j, k))
               * (a(i - 1, j, k) + a(i, j, k));
        result(i, j, k) -= flux / (coords->dx(i, j, k) * coords->J(i, j, k));

        // Flow up
        flux = gU * 0.25 * (coords->J(i, j, k) + coords->J(i, j, kp))
               * (a(i, j, k) + a(i, j, kp));
        result(i, j, k) += flux / (coords->dz(i, j, k) * coords->J(i, j, k));

        // Flow down
        flux = gD * 0.25 * (coords->J(i, j, km) + coords->J(i, j, k))
               * (a(i, j, km) + a(i, j, k));
        result(i, j, k) += flux / (coords->dz(i, j, k) * coords->J(i, j, k));
      }
    }
  }

  return result;
}

} // Namespace FV
