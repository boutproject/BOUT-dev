/**************************************************************************
 * Differential geometry
 * Calculates the covariant metric tensor, and christoffel symbol terms
 * given the contravariant metric tensor terms
 **************************************************************************/

#include <bout/assert.hxx>
#include <bout/constants.hxx>
#include <bout/coordinates.hxx>
#include <msg_stack.hxx>
#include <output.hxx>
#include <utils.hxx>

#include <derivs.hxx>
#include <fft.hxx>
#include <interpolation.hxx>

#include <globals.hxx>

// use anonymous namespace so this utility function is not available outside this file
namespace {
  /// Interpolate a Field2D to a new CELL_LOC with interp_to.
  /// Communicates to set internal guard cells.
  /// Boundary guard cells are set by extrapolating from the grid, like
  /// 'free_o3' boundary conditions
  /// Corner guard cells are set to BoutNaN
  Field2D interpolateAndExtrapolate(const Field2D &f, CELL_LOC location, bool extrap_at_branch_cut) {
    Mesh* localmesh = f.getMesh();
    Field2D result = interp_to(f, location, RGN_NOBNDRY);
    // Ensure result's data is unique. Otherwise result might be a duplicate of
    // f (if no interpolation is needed, e.g. if interpolation is in the
    // z-direction); then f would be communicated. Since this function is used
    // on geometrical quantities that might not be periodic in y even on closed
    // field lines (due to dependence on integrated shear), we don't want to
    // communicate f. We will sort out result's boundary guard cells below, but
    // not f's so we don't want to change f.
    result.allocate();
    localmesh->communicate(result);

    // Extrapolate into boundaries so that differential geometry terms can be
    // interpolated if necessary
    // Note: cannot use applyBoundary("free_o3") here because applyBoundary()
    // would try to create a new Coordinates object since we have not finished
    // initializing yet, leading to an infinite recursion.
    // Also, here we interpolate for the boundary points at xstart/ystart and
    // (xend+1)/(yend+1) instead of extrapolating.
    for (auto bndry : localmesh->getBoundaries()) {
      int extrap_start = 0;
      if ( (location == CELL_XLOW) && (bndry->bx>0) )
        extrap_start = 1;
      else if ( (location == CELL_YLOW) && (bndry->by>0) )
        extrap_start = 1;
      for(bndry->first(); !bndry->isDone(); bndry->next1d()) {
        // interpolate extra boundary point that is missed by interp_to, if
        // necessary
        if (extrap_start>0) {
          // note that either bx or by is >0 here
          result(bndry->x, bndry->y) =
            ( 9.*(f(bndry->x-bndry->bx, bndry->y-bndry->by)
                  + f(bndry->x, bndry->y))
              - f(bndry->x-2*bndry->bx, bndry->y-2*bndry->by)
              - f(bndry->x+bndry->bx, bndry->y+bndry->by)
            )/16.;
        }

        // set boundary guard cells
        if ((bndry->bx != 0 && localmesh->GlobalNx-2*bndry->width >= 3) || (bndry->by != 0 && localmesh->GlobalNy-2*bndry->width >= 3)) {
          if (bndry->bx != 0 && localmesh->LocalNx == 1 && bndry->width == 1) {
            throw BoutException("Not enough points in the x-direction on this "
                "processor for extrapolation needed to use staggered grids. "
                "Increase number of x-guard cells MXG or decrease number of "
                "processors in the x-direction NXPE.");
          }
          if (bndry->by != 0 && localmesh->LocalNy == 1 && bndry->width == 1) {
            throw BoutException("Not enough points in the y-direction on this "
                "processor for extrapolation needed to use staggered grids. "
                "Increase number of y-guard cells MYG or decrease number of "
                "processors in the y-direction NYPE.");
          }
          // extrapolate into boundary guard cells if there are enough grid points
          for(int i=extrap_start;i<bndry->width;i++) {
            int xi = bndry->x + i*bndry->bx;
            int yi = bndry->y + i*bndry->by;
            result(xi, yi) = 3.0*result(xi - bndry->bx, yi - bndry->by) - 3.0*result(xi - 2*bndry->bx, yi - 2*bndry->by) + result(xi - 3*bndry->bx, yi - 3*bndry->by);
          }
        } else {
          // not enough grid points to extrapolate, set equal to last grid point
          for(int i=extrap_start;i<bndry->width;i++) {
            result(bndry->x + i*bndry->bx, bndry->y + i*bndry->by) = result(bndry->x - bndry->bx, bndry->y - bndry->by);
          }
        }
      }
    }

    if (extrap_at_branch_cut) {
      // Extrapolate into guard cells at branch cuts on core/PF field lines
      // This overwrites any communicated values in the guard cells
      int firstjup = location==CELL_YLOW ? localmesh->yend+1 : localmesh->yend;
      for (int i=localmesh->xstart; i<=localmesh->xend; i++) {
        // Lower processor boundary
        if (localmesh->hasBranchCutDown(i)) {
          for (int j=localmesh->ystart-1; j>0; j--) {
            result(i, j) = 3.0*result(i, j+1) - 3.0*result(i, j+2) + result(i, j+3);
          }
        }
        // Upper processor boundary
        if (localmesh->hasBranchCutUp(i)) {
          if (location == CELL_YLOW) {
            // interpolate boundary point, to be symmetric with lower boundary
            int j = localmesh->yend;
            result(i, j) =
              ( 9.*(f(i, j-1) + f(i, j)) - f(i, j-2) - f(i, j+1))/16.;
          }
          for (int j=firstjup; j<localmesh->LocalNy; j++) {
            result(i, j) = 3.0*result(i, j-1) - 3.0*result(i, j-2) + result(i, j-3);
          }
        }
      }
    }

    // Set corner guard cells
    for (int i=0; i<localmesh->xstart; i++) {
      for (int j=0; j<localmesh->ystart; j++) {
        result(i, j) = BoutNaN;
        result(i, localmesh->LocalNy-1-j) = BoutNaN;
        result(localmesh->LocalNx-1-i, j) = BoutNaN;
        result(localmesh->LocalNx-1-i, localmesh->LocalNy-1-j) = BoutNaN;
      }
    }

    return result;
  }
}

Coordinates::Coordinates(Mesh *mesh)
    : dx(1, mesh), dy(1, mesh), dz(1), d1_dx(mesh), d1_dy(mesh), J(1, mesh), Bxy(1, mesh),
      // Identity metric tensor
      g11(1, mesh), g22(1, mesh), g33(1, mesh), g12(0, mesh), g13(0, mesh), g23(0, mesh),
      g_11(1, mesh), g_22(1, mesh), g_33(1, mesh), g_12(0, mesh), g_13(0, mesh),
      g_23(0, mesh), G1_11(mesh), G1_22(mesh), G1_33(mesh), G1_12(mesh), G1_13(mesh),
      G1_23(mesh), G2_11(mesh), G2_22(mesh), G2_33(mesh), G2_12(mesh), G2_13(mesh),
      G2_23(mesh), G3_11(mesh), G3_22(mesh), G3_33(mesh), G3_12(mesh), G3_13(mesh),
      G3_23(mesh), G1(mesh), G2(mesh), G3(mesh), zShift(mesh), ShiftTorsion(mesh),
      IntShiftTorsion(mesh), localmesh(mesh), location(CELL_CENTRE) {

  // If the source does not include y-boundary guard cells then communications
  // are needed to fill them. This provides backwards compatibility even though
  // the communicated guard cell values may not be correct, because the
  // integrated shear is not necessarily continuous when poloidal angle goes
  // from 2pi->0
  bool need_comms = !mesh->sourceHasYGuards();

  if (localmesh->get(dx, "dx", 0., need_comms)) {
    output_warn.write("\tWARNING: differencing quantity 'dx' not found. Set to 1.0\n");
    dx = 1.0;
  }

  if (localmesh->periodicX) {
    localmesh->communicate(dx);
  }

  if (localmesh->get(dy, "dy", 0., need_comms)) {
    output_warn.write("\tWARNING: differencing quantity 'dy' not found. Set to 1.0\n");
    dy = 1.0;
  }

  nz = localmesh->LocalNz;

  if (localmesh->get(dz, "dz")) {
    // Couldn't read dz from input
    int zperiod;
    BoutReal ZMIN, ZMAX;
    Options *options = Options::getRoot();
    if (options->isSet("zperiod")) {
      OPTION(options, zperiod, 1);
      ZMIN = 0.0;
      ZMAX = 1.0 / static_cast<BoutReal>(zperiod);
    } else {
      OPTION(options, ZMIN, 0.0);
      OPTION(options, ZMAX, 1.0);

      zperiod = ROUND(1.0 / (ZMAX - ZMIN));
    }

    dz = (ZMAX - ZMIN) * TWOPI / nz;
  }

  // Diagonal components of metric tensor g^{ij} (default to 1)
  localmesh->get(g11, "g11", 1.0, need_comms);
  localmesh->get(g22, "g22", 1.0, need_comms);
  localmesh->get(g33, "g33", 1.0, need_comms);

  // Off-diagonal elements. Default to 0
  localmesh->get(g12, "g12", 0.0, need_comms);
  localmesh->get(g13, "g13", 0.0, need_comms);
  localmesh->get(g23, "g23", 0.0, need_comms);

  // Check input metrics
  if ((!finite(g11)) || (!finite(g22)) || (!finite(g33))) {
    throw BoutException("\tERROR: Diagonal metrics are not finite!\n");
  }
  if ((min(g11) <= 0.0) || (min(g22) <= 0.0) || (min(g33) <= 0.0)) {
    throw BoutException("\tERROR: Diagonal metrics are negative!\n");
  }
  if ((!finite(g12)) || (!finite(g13)) || (!finite(g23))) {
    throw BoutException("\tERROR: Off-diagonal metrics are not finite!\n");
  }

  /// Find covariant metric components
  // Check if any of the components are present
  if (localmesh->sourceHasVar("g_11") or localmesh->sourceHasVar("g_22") or
      localmesh->sourceHasVar("g_33") or localmesh->sourceHasVar("g_12") or
      localmesh->sourceHasVar("g_13") or localmesh->sourceHasVar("g_23")) {
    // Check that all components are present
    if (localmesh->sourceHasVar("g_11") and localmesh->sourceHasVar("g_22") and
        localmesh->sourceHasVar("g_33") and localmesh->sourceHasVar("g_12") and
        localmesh->sourceHasVar("g_13") and localmesh->sourceHasVar("g_23")) {
      localmesh->get(g_11, "g_11", 0.0, need_comms);
      localmesh->get(g_22, "g_22", 0.0, need_comms);
      localmesh->get(g_33, "g_33", 0.0, need_comms);
      localmesh->get(g_12, "g_12", 0.0, need_comms);
      localmesh->get(g_13, "g_13", 0.0, need_comms);
      localmesh->get(g_23, "g_23", 0.0, need_comms);

      output_warn.write("\tWARNING! Covariant components of metric tensor set manually. "
                        "Contravariant components NOT recalculated\n");

    } else {
      output_warn.write("Not all covariant components of metric tensor found. "
                        "Calculating all from the contravariant tensor\n");
      /// Calculate contravariant metric components if not found
      if (calcCovariant()) {
        throw BoutException("Error in calcCovariant call");
      }
    }
  } else {
    /// Calculate contravariant metric components if not found
    if (calcCovariant()) {
      throw BoutException("Error in calcCovariant call");
    }
  }

  /// Calculate Jacobian and Bxy
  if (jacobian())
    throw BoutException("Error in jacobian call");

  // Attempt to read J from the grid file
  Field2D Jcalc = J;
  if (localmesh->get(J, "J", 0., need_comms)) {
    output_warn.write("\tWARNING: Jacobian 'J' not found. Calculating from metric tensor\n");
    J = Jcalc;
  } else {
    // Compare calculated and loaded values
    output_warn.write("\tMaximum difference in J is %e\n", max(abs(J - Jcalc)));

    // Re-evaluate Bxy using new J
    Bxy = sqrt(g_22) / J;
  }

  // Attempt to read Bxy from the grid file
  Field2D Bcalc = Bxy;
  if (localmesh->get(Bxy, "Bxy", 0., need_comms)) {
    output_warn.write("\tWARNING: Magnitude of B field 'Bxy' not found. Calculating from "
                      "metric tensor\n");
    Bxy = Bcalc;
  } else {
    output_warn.write("\tMaximum difference in Bxy is %e\n", max(abs(Bxy - Bcalc)));
    // Check Bxy
    if (!finite(Bxy)) {
      throw BoutException("\tERROR: Bxy not finite everywhere!\n");
    }
  }

  //////////////////////////////////////////////////////
  /// Calculate Christoffel symbols. Needs communication
  if (geometry()) {
    throw BoutException("Differential geometry failed\n");
  }

  if (localmesh->get(ShiftTorsion, "ShiftTorsion", 0., need_comms)) {
    output_warn.write("\tWARNING: No Torsion specified for zShift. Derivatives may not be correct\n");
    ShiftTorsion = 0.0;
  }

  //////////////////////////////////////////////////////

  // try to read zShift from grid
  if(localmesh->get(zShift, "zShift", 0., need_comms)) {
    // No zShift variable. Try qinty in BOUT grid files
    localmesh->get(zShift, "qinty", 0., need_comms);
  }
  if (need_comms) {
    // zShift should never be periodic on closed field lines, so extrapolate
    // into branch-cut guard cells
    interpolateAndExtrapolate(zShift, CELL_CENTRE, true);
  }

  if (localmesh->IncIntShear) {
    if (localmesh->get(IntShiftTorsion, "IntShiftTorsion", 0., need_comms)) {
      output_warn.write("\tWARNING: No Integrated torsion specified\n");
      IntShiftTorsion = 0.0;
    }
  }
}

Coordinates::Coordinates(Mesh *mesh, const CELL_LOC loc, const Coordinates* coords_in)
    : dx(1, mesh), dy(1, mesh), dz(1), d1_dx(mesh), d1_dy(mesh), J(1, mesh), Bxy(1, mesh),
      // Identity metric tensor
      g11(1, mesh), g22(1, mesh), g33(1, mesh), g12(0, mesh), g13(0, mesh), g23(0, mesh),
      g_11(1, mesh), g_22(1, mesh), g_33(1, mesh), g_12(0, mesh), g_13(0, mesh),
      g_23(0, mesh), G1_11(mesh), G1_22(mesh), G1_33(mesh), G1_12(mesh), G1_13(mesh),
      G1_23(mesh), G2_11(mesh), G2_22(mesh), G2_33(mesh), G2_12(mesh), G2_13(mesh),
      G2_23(mesh), G3_11(mesh), G3_22(mesh), G3_33(mesh), G3_12(mesh), G3_13(mesh),
      G3_23(mesh), G1(mesh), G2(mesh), G3(mesh), zShift(mesh), ShiftTorsion(mesh),
      IntShiftTorsion(mesh), localmesh(mesh), location(loc) {

  dx = interpolateAndExtrapolate(coords_in->dx, location, false);
  dy = interpolateAndExtrapolate(coords_in->dy, location, false);

  nz = localmesh->LocalNz;

  dz = coords_in->dz;

  // Diagonal components of metric tensor g^{ij}
  g11 = interpolateAndExtrapolate(coords_in->g11, location, mesh->hasBranchCut());
  g22 = interpolateAndExtrapolate(coords_in->g22, location, mesh->hasBranchCut());
  g33 = interpolateAndExtrapolate(coords_in->g33, location, mesh->hasBranchCut());

  // Off-diagonal elements.
  g12 = interpolateAndExtrapolate(coords_in->g12, location, mesh->hasBranchCut());
  g13 = interpolateAndExtrapolate(coords_in->g13, location, mesh->hasBranchCut());
  g23 = interpolateAndExtrapolate(coords_in->g23, location, mesh->hasBranchCut());

  // always extrapolate zShift
  zShift = interpolateAndExtrapolate(coords_in->zShift, location, true);

  // Check input metrics
  if ((!finite(g11, RGN_NOBNDRY)) || (!finite(g22, RGN_NOBNDRY)) || (!finite(g33, RGN_NOBNDRY))) {
    throw BoutException("\tERROR: Interpolated diagonal metrics are not finite!\n");
  }
  if ((min(g11) <= 0.0) || (min(g22) <= 0.0) || (min(g33) <= 0.0)) {
    throw BoutException("\tERROR: Interpolated diagonal metrics are negative!\n");
  }
  if ((!finite(g12, RGN_NOBNDRY)) || (!finite(g13, RGN_NOBNDRY)) || (!finite(g23, RGN_NOBNDRY))) {
    throw BoutException("\tERROR: Interpolated off-diagonal metrics are not finite!\n");
  }

  /// Always calculate contravariant metric components so that they are
  /// consistent with the interpolated covariant components
  if (calcCovariant()) {
    throw BoutException("Error in calcCovariant call");
  }

  /// Calculate Jacobian and Bxy
  if (jacobian())
    throw BoutException("Error in jacobian call");

  //////////////////////////////////////////////////////
  /// Calculate Christoffel symbols. Needs communication
  if (geometry()) {
    throw BoutException("Differential geometry failed\n");
  }

  ShiftTorsion = interpolateAndExtrapolate(coords_in->ShiftTorsion, location, false);

  //////////////////////////////////////////////////////

  if (localmesh->IncIntShear) {
    IntShiftTorsion = interpolateAndExtrapolate(coords_in->IntShiftTorsion, location, true);
  }
}

void Coordinates::outputVars(Datafile &file) {
  file.add(dx, "dx", false);
  file.add(dy, "dy", false);
  file.add(dz, "dz", false);

  file.add(g11, "g11", false);
  file.add(g22, "g22", false);
  file.add(g33, "g33", false);
  file.add(g12, "g12", false);
  file.add(g13, "g13", false);
  file.add(g23, "g23", false);

  file.add(g_11, "g_11", false);
  file.add(g_22, "g_22", false);
  file.add(g_33, "g_33", false);
  file.add(g_12, "g_12", false);
  file.add(g_13, "g_13", false);
  file.add(g_23, "g_23", false);

  file.add(J, "J", false);
}

int Coordinates::geometry() {
  TRACE("Coordinates::geometry");

  output_progress.write("Calculating differential geometry terms\n");

  if (min(abs(dx)) < 1e-8)
    throw BoutException("dx magnitude less than 1e-8");

  if (min(abs(dy)) < 1e-8)
    throw BoutException("dy magnitude less than 1e-8");

  if (fabs(dz) < 1e-8)
    throw BoutException("dz magnitude less than 1e-8");

  // Check input metrics
  if ((!finite(g11, RGN_NOBNDRY)) || (!finite(g22, RGN_NOBNDRY)) || (!finite(g33, RGN_NOBNDRY))) {
    throw BoutException("\tERROR: Diagonal metrics are not finite!\n");
  }
  if ((min(g11) <= 0.0) || (min(g22) <= 0.0) || (min(g33) <= 0.0)) {
    throw BoutException("\tERROR: Diagonal metrics are negative!\n");
  }
  if ((!finite(g12, RGN_NOBNDRY)) || (!finite(g13, RGN_NOBNDRY)) || (!finite(g23, RGN_NOBNDRY))) {
    throw BoutException("\tERROR: Off-diagonal metrics are not finite!\n");
  }

  if ((!finite(g_11, RGN_NOBNDRY)) || (!finite(g_22, RGN_NOBNDRY)) || (!finite(g_33, RGN_NOBNDRY))) {
    throw BoutException("\tERROR: Diagonal g_ij metrics are not finite!\n");
  }
  if ((min(g_11) <= 0.0) || (min(g_22) <= 0.0) || (min(g_33) <= 0.0)) {
    throw BoutException("\tERROR: Diagonal g_ij metrics are negative!\n");
  }
  if ((!finite(g_12, RGN_NOBNDRY)) || (!finite(g_13, RGN_NOBNDRY)) || (!finite(g_23, RGN_NOBNDRY))) {
    throw BoutException("\tERROR: Off-diagonal g_ij metrics are not finite!\n");
  }

  // Calculate Christoffel symbol terms (18 independent values)
  // Note: This calculation is completely general: metric
  // tensor can be 2D or 3D. For 2D, all DDZ terms are zero

  G1_11 = 0.5 * g11 * DDX(g_11) + g12 * (DDX(g_12) - 0.5 * DDY(g_11)) +
          g13 * (DDX(g_13) - 0.5 * DDZ(g_11));
  G1_22 = g11 * (DDY(g_12) - 0.5 * DDX(g_22)) + 0.5 * g12 * DDY(g_22) +
          g13 * (DDY(g_23) - 0.5 * DDZ(g_22));
  G1_33 = g11 * (DDZ(g_13) - 0.5 * DDX(g_33)) + g12 * (DDZ(g_23) - 0.5 * DDY(g_33)) +
          0.5 * g13 * DDZ(g_33);
  G1_12 = 0.5 * g11 * DDY(g_11) + 0.5 * g12 * DDX(g_22) +
          0.5 * g13 * (DDY(g_13) + DDX(g_23) - DDZ(g_12));
  G1_13 = 0.5 * g11 * DDZ(g_11) + 0.5 * g12 * (DDZ(g_12) + DDX(g_23) - DDY(g_13)) +
          0.5 * g13 * DDX(g_33);
  G1_23 = 0.5 * g11 * (DDZ(g_12) + DDY(g_13) - DDX(g_23)) +
          0.5 * g12 * (DDZ(g_22) + DDY(g_23) - DDY(g_23))
          // + 0.5 *g13*(DDZ(g_32) + DDY(g_33) - DDZ(g_23));
          // which equals
          + 0.5 * g13 * DDY(g_33);

  G2_11 = 0.5 * g12 * DDX(g_11) + g22 * (DDX(g_12) - 0.5 * DDY(g_11)) +
          g23 * (DDX(g_13) - 0.5 * DDZ(g_11));
  G2_22 = g12 * (DDY(g_12) - 0.5 * DDX(g_22)) + 0.5 * g22 * DDY(g_22) +
          g23 * (DDY(g23) - 0.5 * DDZ(g_22));
  G2_33 = g12 * (DDZ(g_13) - 0.5 * DDX(g_33)) + g22 * (DDZ(g_23) - 0.5 * DDY(g_33)) +
          0.5 * g23 * DDZ(g_33);
  G2_12 = 0.5 * g12 * DDY(g_11) + 0.5 * g22 * DDX(g_22) +
          0.5 * g23 * (DDY(g_13) + DDX(g_23) - DDZ(g_12));
  G2_13 =
      // 0.5 *g21*(DDZ(g_11) + DDX(g_13) - DDX(g_13))
      // which equals
      0.5 * g12 * (DDZ(g_11) + DDX(g_13) - DDX(g_13))
      // + 0.5 *g22*(DDZ(g_21) + DDX(g_23) - DDY(g_13))
      // which equals
      + 0.5 * g22 * (DDZ(g_12) + DDX(g_23) - DDY(g_13))
      // + 0.5 *g23*(DDZ(g_31) + DDX(g_33) - DDZ(g_13));
      // which equals
      + 0.5 * g23 * DDX(g_33);
  G2_23 = 0.5 * g12 * (DDZ(g_12) + DDY(g_13) - DDX(g_23)) + 0.5 * g22 * DDZ(g_22) +
          0.5 * g23 * DDY(g_33);

  G3_11 = 0.5 * g13 * DDX(g_11) + g23 * (DDX(g_12) - 0.5 * DDY(g_11)) +
          g33 * (DDX(g_13) - 0.5 * DDZ(g_11));
  G3_22 = g13 * (DDY(g_12) - 0.5 * DDX(g_22)) + 0.5 * g23 * DDY(g_22) +
          g33 * (DDY(g_23) - 0.5 * DDZ(g_22));
  G3_33 = g13 * (DDZ(g_13) - 0.5 * DDX(g_33)) + g23 * (DDZ(g_23) - 0.5 * DDY(g_33)) +
          0.5 * g33 * DDZ(g_33);
  G3_12 =
      // 0.5 *g31*(DDY(g_11) + DDX(g_12) - DDX(g_12))
      // which equals to
      0.5 * g13 * DDY(g_11)
      // + 0.5 *g32*(DDY(g_21) + DDX(g_22) - DDY(g_12))
      // which equals to
      + 0.5 * g23 * DDX(g_22)
      //+ 0.5 *g33*(DDY(g_31) + DDX(g_32) - DDZ(g_12));
      // which equals to
      + 0.5 * g33 * (DDY(g_13) + DDX(g_23) - DDZ(g_12));
  G3_13 = 0.5 * g13 * DDZ(g_11) + 0.5 * g23 * (DDZ(g_12) + DDX(g_23) - DDY(g_13)) +
          0.5 * g33 * DDX(g_33);
  G3_23 = 0.5 * g13 * (DDZ(g_12) + DDY(g_13) - DDX(g_23)) + 0.5 * g23 * DDZ(g_22) +
          0.5 * g33 * DDY(g_33);

  G1 = (DDX(J * g11) + DDY(J * g12) + DDZ(J * g13)) / J;
  G2 = (DDX(J * g12) + DDY(J * g22) + DDZ(J * g23)) / J;
  G3 = (DDX(J * g13) + DDY(J * g23) + DDZ(J * g33)) / J;

  // Communicate christoffel symbol terms
  output_progress.write("\tCommunicating connection terms\n");

  FieldGroup com;

  com.add(G1_11);
  com.add(G1_22);
  com.add(G1_33);
  com.add(G1_12);
  com.add(G1_13);
  com.add(G1_23);

  com.add(G2_11);
  com.add(G2_22);
  com.add(G2_33);
  com.add(G2_12);
  com.add(G2_13);
  com.add(G2_23);

  com.add(G3_11);
  com.add(G3_22);
  com.add(G3_33);
  com.add(G3_12);
  com.add(G3_13);
  com.add(G3_23);

  com.add(G1);
  com.add(G2);
  com.add(G3);

  localmesh->communicate(com);

  //////////////////////////////////////////////////////
  /// Non-uniform meshes. Need to use DDX, DDY

  OPTION(Options::getRoot(), non_uniform, true);

  // If the source does not include y-boundary guard cells then communications
  // are needed to fill them. This provides backwards compatibility even though
  // the communicated guard cell values may not be correct, because the
  // integrated shear is not necessarily continuous when poloidal angle goes
  // from 2pi->0
  bool need_comms = !mesh->sourceHasYGuards();

  Field2D d2x, d2y; // d^2 x / d i^2
  // Read correction for non-uniform meshes
  if (localmesh->get(d2x, "d2x", 0., need_comms)) {
    output_warn.write(
        "\tWARNING: differencing quantity 'd2x' not found. Calculating from dx\n");
    d1_dx = localmesh->indexDDX(1. / dx); // d/di(1/dx)
  } else {
    d1_dx = -d2x / (dx * dx);
  }

  if (localmesh->get(d2y, "d2y", 0., need_comms)) {
    output_warn.write(
        "\tWARNING: differencing quantity 'd2y' not found. Calculating from dy\n");
    d1_dy = localmesh->indexDDY(1. / dy); // d/di(1/dy)
  } else {
    d1_dy = -d2y / (dy * dy);
  }

  return 0;
}

int Coordinates::calcCovariant() {
  TRACE("Coordinates::calcCovariant");

  // Make sure metric elements are allocated
  g_11.allocate();
  g_22.allocate();
  g_33.allocate();
  g_12.allocate();
  g_13.allocate();
  g_23.allocate();

  g_11.setLocation(location);
  g_22.setLocation(location);
  g_33.setLocation(location);
  g_12.setLocation(location);
  g_13.setLocation(location);
  g_23.setLocation(location);

  // Perform inversion of g^{ij} to get g_{ij}
  // NOTE: Currently this bit assumes that metric terms are Field2D objects

  auto a = Matrix<BoutReal>(3, 3);

  for (int jx = 0; jx < localmesh->LocalNx; jx++) {
    for (int jy = 0; jy < localmesh->LocalNy; jy++) {
      // set elements of g
      a(0, 0) = g11(jx, jy);
      a(1, 1) = g22(jx, jy);
      a(2, 2) = g33(jx, jy);

      a(0, 1) = a(1, 0) = g12(jx, jy);
      a(1, 2) = a(2, 1) = g23(jx, jy);
      a(0, 2) = a(2, 0) = g13(jx, jy);

      // invert
      if (invert3x3(a)) {
        output_error.write("\tERROR: metric tensor is singular at (%d, %d)\n", jx, jy);
        return 1;
      }

      // put elements into g_{ij}
      g_11(jx, jy) = a(0, 0);
      g_22(jx, jy) = a(1, 1);
      g_33(jx, jy) = a(2, 2);

      g_12(jx, jy) = a(0, 1);
      g_13(jx, jy) = a(0, 2);
      g_23(jx, jy) = a(1, 2);
    }
  }

  BoutReal maxerr;
  maxerr = BOUTMAX(max(abs((g_11 * g11 + g_12 * g12 + g_13 * g13) - 1)),
                   max(abs((g_12 * g12 + g_22 * g22 + g_23 * g23) - 1)),
                   max(abs((g_13 * g13 + g_23 * g23 + g_33 * g33) - 1)));

  output_info.write("\tLocal maximum error in diagonal inversion is %e\n", maxerr);

  maxerr = BOUTMAX(max(abs(g_11 * g12 + g_12 * g22 + g_13 * g23)),
                   max(abs(g_11 * g13 + g_12 * g23 + g_13 * g33)),
                   max(abs(g_12 * g13 + g_22 * g23 + g_23 * g33)));

  output_info.write("\tLocal maximum error in off-diagonal inversion is %e\n", maxerr);

  return 0;
}

int Coordinates::calcContravariant() {
  TRACE("Coordinates::calcContravariant");

  // Make sure metric elements are allocated
  g11.allocate();
  g22.allocate();
  g33.allocate();
  g12.allocate();
  g13.allocate();
  g23.allocate();

  // Perform inversion of g_{ij} to get g^{ij}
  // NOTE: Currently this bit assumes that metric terms are Field2D objects

  auto a = Matrix<BoutReal>(3, 3);

  for (int jx = 0; jx < localmesh->LocalNx; jx++) {
    for (int jy = 0; jy < localmesh->LocalNy; jy++) {
      // set elements of g
      a(0, 0) = g_11(jx, jy);
      a(1, 1) = g_22(jx, jy);
      a(2, 2) = g_33(jx, jy);

      a(0, 1) = a(1, 0) = g_12(jx, jy);
      a(1, 2) = a(2, 1) = g_23(jx, jy);
      a(0, 2) = a(2, 0) = g_13(jx, jy);

      // invert
      if (invert3x3(a)) {
        output_error.write("\tERROR: metric tensor is singular at (%d, %d)\n", jx, jy);
        return 1;
      }

      // put elements into g_{ij}
      g11(jx, jy) = a(0, 0);
      g22(jx, jy) = a(1, 1);
      g33(jx, jy) = a(2, 2);

      g12(jx, jy) = a(0, 1);
      g13(jx, jy) = a(0, 2);
      g23(jx, jy) = a(1, 2);
    }
  }

  BoutReal maxerr;
  maxerr = BOUTMAX(max(abs((g_11 * g11 + g_12 * g12 + g_13 * g13) - 1)),
                   max(abs((g_12 * g12 + g_22 * g22 + g_23 * g23) - 1)),
                   max(abs((g_13 * g13 + g_23 * g23 + g_33 * g33) - 1)));

  output_info.write("\tMaximum error in diagonal inversion is %e\n", maxerr);

  maxerr = BOUTMAX(max(abs(g_11 * g12 + g_12 * g22 + g_13 * g23)),
                   max(abs(g_11 * g13 + g_12 * g23 + g_13 * g33)),
                   max(abs(g_12 * g13 + g_22 * g23 + g_23 * g33)));

  output_info.write("\tMaximum error in off-diagonal inversion is %e\n", maxerr);
  return 0;
}

int Coordinates::jacobian() {
  TRACE("Coordinates::jacobian");
  // calculate Jacobian using g^-1 = det[g^ij], J = sqrt(g)

  Field2D g = g11 * g22 * g33 + 2.0 * g12 * g13 * g23 - g11 * g23 * g23 -
              g22 * g13 * g13 - g33 * g12 * g12;

  // Check that g is positive
  if (min(g) < 0.0) {
    throw BoutException("The determinant of g^ij is somewhere less than 0.0");
  }
  J = 1. / sqrt(g);

  // Check jacobian
  if (!finite(J, RGN_NOBNDRY)) {
    throw BoutException("\tERROR: Jacobian not finite everywhere!\n");
  }
  if (min(abs(J)) < 1.0e-10) {
    throw BoutException("\tERROR: Jacobian becomes very small\n");
  }

  if (min(g_22) < 0.0) {
    throw BoutException("g_22 is somewhere less than 0.0");
  }
  Bxy = sqrt(g_22) / J;

  return 0;
}

/*******************************************************************************
 * Operators
 *
 *******************************************************************************/

const Field2D Coordinates::DDX(const Field2D &f, CELL_LOC loc, DIFF_METHOD method, REGION region) {
  ASSERT1(location == loc || loc == CELL_DEFAULT);
  return localmesh->indexDDX(f, loc, method, region) / dx;
}

const Field2D Coordinates::DDY(const Field2D &f, CELL_LOC loc, DIFF_METHOD method, REGION region) {
  ASSERT1(location == loc || loc == CELL_DEFAULT);
  return localmesh->indexDDY(f, loc, method, region) / dy;
}

const Field2D Coordinates::DDZ(const Field2D &f, CELL_LOC loc,
                               DIFF_METHOD UNUSED(method), REGION UNUSED(region)) {
  ASSERT1(location == loc || loc == CELL_DEFAULT);
  ASSERT1(f.getMesh() == localmesh);
  auto result = Field2D(0.0, localmesh);
  result.setLocation(location);
  return result;
}

#include <derivs.hxx>

/////////////////////////////////////////////////////////
// Parallel gradient

const Field2D Coordinates::Grad_par(const Field2D &var, CELL_LOC outloc,
                                    DIFF_METHOD UNUSED(method)) {
  TRACE("Coordinates::Grad_par( Field2D )");
  ASSERT1(location == outloc || outloc == CELL_DEFAULT);

  return DDY(var) / sqrt(g_22);
}

const Field3D Coordinates::Grad_par(const Field3D &var, CELL_LOC outloc,
                                    DIFF_METHOD method) {
  TRACE("Coordinates::Grad_par( Field3D )");
  ASSERT1(location == outloc || outloc == CELL_DEFAULT);

  return ::DDY(var, outloc, method) / sqrt(g_22);
}

/////////////////////////////////////////////////////////
// Vpar_Grad_par
// vparallel times the parallel derivative along unperturbed B-field

const Field2D Coordinates::Vpar_Grad_par(const Field2D &v, const Field2D &f,
                                         CELL_LOC outloc,
                                         DIFF_METHOD UNUSED(method)) {
  ASSERT1(location == outloc || outloc == CELL_DEFAULT);
  return VDDY(v, f) / sqrt(g_22);
}

const Field3D Coordinates::Vpar_Grad_par(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                                         DIFF_METHOD method) {
  ASSERT1(location == outloc || outloc == CELL_DEFAULT);
  return VDDY(v, f, outloc, method) / sqrt(g_22);
}

/////////////////////////////////////////////////////////
// Parallel divergence

const Field2D Coordinates::Div_par(const Field2D &f, CELL_LOC outloc,
                                   DIFF_METHOD method) {
  TRACE("Coordinates::Div_par( Field2D )");
  ASSERT1(location == outloc || outloc == CELL_DEFAULT);

  // Need Bxy at location of f, which might be different from location of this
  // Coordinates object
  Field2D Bxy_floc = f.getCoordinates()->Bxy;

  return Bxy * Grad_par(f / Bxy_floc, outloc, method);
}

const Field3D Coordinates::Div_par(const Field3D &f, CELL_LOC outloc,
                                   DIFF_METHOD method) {
  TRACE("Coordinates::Div_par( Field3D )");
  ASSERT1(location == outloc || outloc == CELL_DEFAULT);
  
  // Need Bxy at location of f, which might be different from location of this
  // Coordinates object
  Field2D Bxy_floc = f.getCoordinates()->Bxy;

  if (f.hasYupYdown()) {
    // Need to modify yup and ydown fields
    Field3D f_B = f / Bxy_floc;
    if (&f.yup() == &f) {
      // Identity, yup and ydown point to same field
      f_B.mergeYupYdown();
    } else {
      // Distinct fields
      f_B.splitYupYdown();
      f_B.yup() = f.yup() / Bxy_floc;
      f_B.ydown() = f.ydown() / Bxy_floc;
    }
    return Bxy * Grad_par(f_B, outloc, method);
  }

  // No yup/ydown fields. The Grad_par operator will
  // shift to field aligned coordinates
  return Bxy * Grad_par(f / Bxy_floc, outloc, method);
}

/////////////////////////////////////////////////////////
// second parallel derivative (b dot Grad)(b dot Grad)
// Note: For parallel Laplacian use Laplace_par

const Field2D Coordinates::Grad2_par2(const Field2D &f, CELL_LOC outloc) {
  TRACE("Coordinates::Grad2_par2( Field2D )");
  ASSERT1(location == outloc || outloc == CELL_DEFAULT);

  Field2D sg = sqrt(g_22);
  Field2D result = DDY(1. / sg, outloc) * DDY(f, outloc) / sg + D2DY2(f, outloc) / g_22;

  return result;
}

const Field3D Coordinates::Grad2_par2(const Field3D &f, CELL_LOC outloc) {
  TRACE("Coordinates::Grad2_par2( Field3D )");
  ASSERT1(location == outloc || outloc == CELL_DEFAULT);

  Field2D sg(localmesh);
  Field3D result(localmesh), r2(localmesh);

  sg = sqrt(g_22);
  sg = DDY(1. / sg) / sg;

  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }

  if (sg.getLocation() != outloc) {
    localmesh->communicate(sg);
    sg = interp_to(sg, outloc);
  }

  result = ::DDY(f, outloc);

  r2 = D2DY2(f, outloc) / interp_to(g_22, outloc);

  result = sg * result + r2;

  ASSERT2(((outloc == CELL_DEFAULT) && (result.getLocation() == f.getLocation())) ||
          (result.getLocation() == outloc));

  return result;
}

/////////////////////////////////////////////////////////
// perpendicular Laplacian operator

#include <invert_laplace.hxx> // Delp2 uses same coefficients as inversion code

const Field2D Coordinates::Delp2(const Field2D &f, CELL_LOC outloc) {
  TRACE("Coordinates::Delp2( Field2D )");
  ASSERT1(location == outloc || outloc == CELL_DEFAULT);

  Field2D result = G1 * DDX(f, outloc) + g11 * D2DX2(f, outloc);

  return result;
}

const Field3D Coordinates::Delp2(const Field3D &f, CELL_LOC outloc) {
  TRACE("Coordinates::Delp2( Field3D )");
  ASSERT1(location == outloc || outloc == CELL_DEFAULT);

  if (localmesh->GlobalNx == 1 && localmesh->GlobalNz == 1) {
    // copy mesh, location, etc
    return f*0;
  }
  ASSERT2(localmesh->xstart > 0); // Need at least one guard cell

  if (outloc == CELL_DEFAULT) outloc = f.getLocation();
  ASSERT2(f.getLocation() == outloc);

  Field3D result(localmesh);
  result.allocate();
  result.setLocation(f.getLocation());

  int ncz = localmesh->LocalNz;

  // Allocate memory
  auto ft = Matrix<dcomplex>(localmesh->LocalNx, ncz / 2 + 1);
  auto delft = Matrix<dcomplex>(localmesh->LocalNx, ncz / 2 + 1);

  // Loop over all y indices
  for (int jy = 0; jy < localmesh->LocalNy; jy++) {

    // Take forward FFT

    for (int jx = 0; jx < localmesh->LocalNx; jx++)
      rfft(&f(jx, jy, 0), ncz, &ft(jx, 0));

    // Loop over kz
    for (int jz = 0; jz <= ncz / 2; jz++) {
      dcomplex a, b, c;

      // No smoothing in the x direction
      for (int jx = localmesh->xstart; jx <= localmesh->xend; jx++) {
        // Perform x derivative

        laplace_tridag_coefs(jx, jy, jz, a, b, c, nullptr, nullptr, outloc);

        delft(jx, jz) = a * ft(jx - 1, jz) + b * ft(jx, jz) + c * ft(jx + 1, jz);
      }
    }

    // Reverse FFT
    for (int jx = localmesh->xstart; jx <= localmesh->xend; jx++) {

      irfft(&delft(jx, 0), ncz, &result(jx, jy, 0));
    }

    // Boundaries
    for (int jz = 0; jz < ncz; jz++) {
      for (int jx = 0; jx < localmesh->xstart; jx++) {
        result(jx, jy, jz) = 0.0;
      }
      for (int jx = localmesh->xend + 1; jx < localmesh->LocalNx; jx++) {
        result(jx, jy, jz) = 0.0;
      }
    }
  }

  ASSERT2(result.getLocation() == f.getLocation());

  return result;
}

const FieldPerp Coordinates::Delp2(const FieldPerp &f, CELL_LOC outloc) {
  TRACE("Coordinates::Delp2( FieldPerp )");

  if (outloc == CELL_DEFAULT) outloc = f.getLocation();

  ASSERT1(location == outloc || outloc == CELL_DEFAULT);
  ASSERT2(f.getLocation() == outloc);

  FieldPerp result(localmesh);
  result.allocate();
  result.setLocation(outloc);

  int jy = f.getIndex();
  result.setIndex(jy);

  int ncz = localmesh->LocalNz;

  // Allocate memory
  auto ft = Matrix<dcomplex>(localmesh->LocalNx, ncz / 2 + 1);
  auto delft = Matrix<dcomplex>(localmesh->LocalNx, ncz / 2 + 1);

  // Take forward FFT
  for (int jx = 0; jx < localmesh->LocalNx; jx++)
    rfft(&f(jx, 0), ncz, &ft(jx, 0));

  // Loop over kz
  for (int jz = 0; jz <= ncz / 2; jz++) {

    // No smoothing in the x direction
    for (int jx = 2; jx < (localmesh->LocalNx - 2); jx++) {
      // Perform x derivative

      dcomplex a, b, c;
      laplace_tridag_coefs(jx, jy, jz, a, b, c);

      delft(jx, jz) = a * ft(jx - 1, jz) + b * ft(jx, jz) + c * ft(jx + 1, jz);
    }
  }

  // Reverse FFT
  for (int jx = 1; jx < (localmesh->LocalNx - 1); jx++) {
    irfft(&delft(jx, 0), ncz, &result(jx, 0));
  }

  // Boundaries
  for (int jz = 0; jz < ncz; jz++) {
    result(0, jz) = 0.0;
    result(localmesh->LocalNx - 1, jz) = 0.0;
  }

  return result;
}

const Field2D Coordinates::Laplace_par(const Field2D &f, CELL_LOC outloc) {
  ASSERT1(location == outloc || outloc == CELL_DEFAULT);
  return D2DY2(f, outloc) / g_22 + DDY(J / g_22, outloc) * DDY(f, outloc) / J;
}

const Field3D Coordinates::Laplace_par(const Field3D &f, CELL_LOC outloc) {
  ASSERT1(location == outloc || outloc == CELL_DEFAULT);
  return D2DY2(f, outloc) / g_22 + DDY(J / g_22, outloc) * ::DDY(f, outloc) / J;
}

// Full Laplacian operator on scalar field

const Field2D Coordinates::Laplace(const Field2D &f, CELL_LOC outloc) {
  TRACE("Coordinates::Laplace( Field2D )");
  ASSERT1(location == outloc || outloc == CELL_DEFAULT);

  Field2D result =
      G1 * DDX(f, outloc) + G2 * DDY(f, outloc) + g11 * D2DX2(f, outloc) + g22 * D2DY2(f, outloc) + 2.0 * g12 * D2DXDY(f, outloc);

  return result;
}

const Field3D Coordinates::Laplace(const Field3D &f, CELL_LOC outloc) {
  TRACE("Coordinates::Laplace( Field3D )");
  ASSERT1(location == outloc || outloc == CELL_DEFAULT);

  Field3D result = G1 * ::DDX(f, outloc) + G2 * ::DDY(f, outloc) + G3 * ::DDZ(f, outloc) + g11 * D2DX2(f, outloc) +
                   g22 * D2DY2(f, outloc) + g33 * D2DZ2(f, outloc) +
                   2.0 * (g12 * D2DXDY(f, outloc) + g13 * D2DXDZ(f, outloc) + g23 * D2DYDZ(f, outloc));

  ASSERT2(result.getLocation() == f.getLocation());

  return result;
}
