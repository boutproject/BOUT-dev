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

#include "parallel/fci.hxx"

// use anonymous namespace so this utility function is not available outside this file
namespace {
/// Interpolate a Field2D to a new CELL_LOC with interp_to.
/// Communicates to set internal guard cells.
/// Boundary guard cells are set by extrapolating from the grid, like
/// 'free_o3' boundary conditions
/// Corner guard cells are set to BoutNaN
Field2D interpolateAndExtrapolate(const Field2D& f, CELL_LOC location,
    bool extrapolate_x = true, bool extrapolate_y = true,
    bool no_extra_interpolate = false) {

  Mesh* localmesh = f.getMesh();
  Field2D result = interp_to(f, location, "RGN_NOBNDRY");
  // Ensure result's data is unique. Otherwise result might be a duplicate of
  // f (if no interpolation is needed, e.g. if interpolation is in the
  // z-direction); then f would be communicated. Since this function is used
  // on geometrical quantities that might not be periodic in y even on closed
  // field lines (due to dependence on integrated shear), we don't want to
  // communicate f. We will sort out result's boundary guard cells below, but
  // not f's so we don't want to change f.
  result.allocate();
  localmesh->communicate(result);

  // Extrapolate into boundaries (if requested) so that differential geometry
  // terms can be interpolated if necessary
  // Note: cannot use applyBoundary("free_o3") here because applyBoundary()
  // would try to create a new Coordinates object since we have not finished
  // initializing yet, leading to an infinite recursion.
  // Also, here we interpolate for the boundary points at xstart/ystart and
  // (xend+1)/(yend+1) instead of extrapolating.
  for (auto& bndry : localmesh->getBoundaries()) {
    if ((extrapolate_x and bndry->bx != 0) or (extrapolate_y and bndry->by != 0)) {
      int extrap_start = 0;
      if (not no_extra_interpolate) {
        // Can use no_extra_interpolate argument to skip the extra interpolation when we
        // want to extrapolate the Christoffel symbol terms which come from derivatives so
        // don't have the extra point set already
        if ((location == CELL_XLOW) && (bndry->bx > 0)) {
          extrap_start = 1;
        } else if ((location == CELL_YLOW) && (bndry->by > 0)) {
          extrap_start = 1;
        }
      }
      for (bndry->first(); !bndry->isDone(); bndry->next1d()) {
        // interpolate extra boundary point that is missed by interp_to, if
        // necessary.
        // Only interpolate this point if we are actually changing location. E.g.
        // when we use this function to extrapolate J and Bxy on staggered grids,
        // this point should already be set correctly because the metric
        // components have been interpolated to here.
        if (extrap_start > 0 and f.getLocation() != location) {
          ASSERT1(bndry->bx == 0 or localmesh->xstart > 1);
          ASSERT1(bndry->by == 0 or localmesh->ystart > 1);
          // note that either bx or by is >0 here
          result(bndry->x, bndry->y) =
              (9. * (f(bndry->x - bndry->bx, bndry->y - bndry->by) + f(bndry->x, bndry->y))
               - f(bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by)
               - f(bndry->x + bndry->bx, bndry->y + bndry->by))
              / 16.;
        }

        // set boundary guard cells
        if ((bndry->bx != 0 && localmesh->GlobalNx - 2 * bndry->width >= 3)
            || (bndry->by != 0
                && localmesh->GlobalNy - localmesh->numberOfYBoundaries() * bndry->width
                   >= 3))
        {
          if (bndry->bx != 0 && localmesh->LocalNx == 1 && bndry->width == 1) {
            throw BoutException(
                "Not enough points in the x-direction on this "
                "processor for extrapolation needed to use staggered grids. "
                "Increase number of x-guard cells MXG or decrease number of "
                "processors in the x-direction NXPE.");
          }
          if (bndry->by != 0 && localmesh->LocalNy == 1 && bndry->width == 1) {
            throw BoutException(
                "Not enough points in the y-direction on this "
                "processor for extrapolation needed to use staggered grids. "
                "Increase number of y-guard cells MYG or decrease number of "
                "processors in the y-direction NYPE.");
          }
          // extrapolate into boundary guard cells if there are enough grid points
          for (int i = extrap_start; i < bndry->width; i++) {
            int xi = bndry->x + i * bndry->bx;
            int yi = bndry->y + i * bndry->by;
            result(xi, yi) = 3.0 * result(xi - bndry->bx, yi - bndry->by)
                             - 3.0 * result(xi - 2 * bndry->bx, yi - 2 * bndry->by)
                             + result(xi - 3 * bndry->bx, yi - 3 * bndry->by);
          }
        } else {
          // not enough grid points to extrapolate, set equal to last grid point
          for (int i = extrap_start; i < bndry->width; i++) {
            result(bndry->x + i * bndry->bx, bndry->y + i * bndry->by) =
                result(bndry->x - bndry->bx, bndry->y - bndry->by);
          }
        }
      }
    }
  }

  // Set corner guard cells
  for (int i = 0; i < localmesh->xstart; i++) {
    for (int j = 0; j < localmesh->ystart; j++) {
      result(i, j) = BoutNaN;
      result(i, localmesh->LocalNy - 1 - j) = BoutNaN;
      result(localmesh->LocalNx - 1 - i, j) = BoutNaN;
      result(localmesh->LocalNx - 1 - i, localmesh->LocalNy - 1 - j) = BoutNaN;
    }
  }

  return result;
}

// If the CELL_CENTRE variable was read, the staggered version is required to
// also exist for consistency
void checkStaggeredGet(Mesh* mesh, const std::string& name, const std::string& suffix) {
  if (mesh->sourceHasVar(name) != mesh->sourceHasVar(name+suffix)) {
    throw BoutException("Attempting to read staggered fields from grid, but " + name
        + " is not present in both CELL_CENTRE and staggered versions.");
  }
}

// convenience function for repeated code
void getAtLoc(Mesh* mesh, Field2D &var, const std::string& name,
    const std::string& suffix, CELL_LOC location, BoutReal default_value = 0.) {

  checkStaggeredGet(mesh, name, suffix);
  mesh->get(var, name+suffix, default_value);
  var.setLocation(location);
}

std::string getLocationSuffix(CELL_LOC location) {
  switch (location) {
  case CELL_CENTRE: {
      return "";
    }
  case CELL_XLOW: {
      return "_xlow";
    }
  case CELL_YLOW: {
      return "_ylow";
    }
  case CELL_ZLOW: {
      // geometrical quantities are Field2D, so CELL_ZLOW version is the same
      // as CELL_CENTRE
      return "";
    }
  default: {
      throw BoutException("Incorrect location passed to "
          "Coordinates(Mesh*,const CELL_LOC,const Coordinates*) constructor.");
    }
  }
}
}

Coordinates::Coordinates(Mesh* mesh, Field2D dx, Field2D dy, BoutReal dz, Field2D J,
                         Field2D Bxy, Field2D g11, Field2D g22, Field2D g33, Field2D g12,
                         Field2D g13, Field2D g23, Field2D g_11, Field2D g_22,
                         Field2D g_33, Field2D g_12, Field2D g_13, Field2D g_23,
                         Field2D ShiftTorsion, Field2D IntShiftTorsion,
                         bool calculate_geometry)
    : dx(std::move(dx)), dy(std::move(dy)), dz(dz), J(std::move(J)), Bxy(std::move(Bxy)),
      g11(std::move(g11)), g22(std::move(g22)), g33(std::move(g33)), g12(std::move(g12)),
      g13(std::move(g13)), g23(std::move(g23)), g_11(std::move(g_11)),
      g_22(std::move(g_22)), g_33(std::move(g_33)), g_12(std::move(g_12)),
      g_13(std::move(g_13)), g_23(std::move(g_23)), ShiftTorsion(std::move(ShiftTorsion)),
      IntShiftTorsion(std::move(IntShiftTorsion)), nz(mesh->LocalNz), localmesh(mesh),
      location(CELL_CENTRE) {
  if (calculate_geometry) {
    if (geometry()) {
      throw BoutException("Differential geometry failed\n");
    }
  }
}

Coordinates::Coordinates(Mesh* mesh, Options* options)
    : dx(1, mesh), dy(1, mesh), dz(1), d1_dx(mesh), d1_dy(mesh), J(1, mesh), Bxy(1, mesh),
      // Identity metric tensor
      g11(1, mesh), g22(1, mesh), g33(1, mesh), g12(0, mesh), g13(0, mesh), g23(0, mesh),
      g_11(1, mesh), g_22(1, mesh), g_33(1, mesh), g_12(0, mesh), g_13(0, mesh),
      g_23(0, mesh), G1_11(mesh), G1_22(mesh), G1_33(mesh), G1_12(mesh), G1_13(mesh),
      G1_23(mesh), G2_11(mesh), G2_22(mesh), G2_33(mesh), G2_12(mesh), G2_13(mesh),
      G2_23(mesh), G3_11(mesh), G3_22(mesh), G3_33(mesh), G3_12(mesh), G3_13(mesh),
      G3_23(mesh), G1(mesh), G2(mesh), G3(mesh), ShiftTorsion(mesh),
      IntShiftTorsion(mesh), localmesh(mesh), location(CELL_CENTRE) {

  if (options == nullptr) {
    options = Options::getRoot()->getSection("mesh");
  }

  // Note: If boundary cells were not loaded from the grid file, use
  // 'interpolateAndExtrapolate' to set them. Ensures that derivatives are
  // smooth at all the boundaries.

  const bool extrapolate_x = (*options)["extrapolate_x"].withDefault(not mesh->sourceHasXBoundaryGuards());
  const bool extrapolate_y = (*options)["extrapolate_y"].withDefault(not mesh->sourceHasYBoundaryGuards());

  if (extrapolate_x) {
    output_warn.write(_("WARNING: extrapolating input mesh quantities into x-boundary "
          "cells. Set option extrapolate_x=false to disable this.\n"));
  }

  if (extrapolate_y) {
    output_warn.write(_("WARNING: extrapolating input mesh quantities into y-boundary "
          "cells. Set option extrapolate_y=false to disable this.\n"));
  }

  mesh->get(dx, "dx", 1.0);
  dx = interpolateAndExtrapolate(dx, location, extrapolate_x, extrapolate_y);

  if (mesh->periodicX) {
    mesh->communicate(dx);
  }

  mesh->get(dy, "dy", 1.0);
  dy = interpolateAndExtrapolate(dy, location, extrapolate_x, extrapolate_y);

  nz = mesh->LocalNz;

  {
    auto& options = Options::root();
    const bool has_zperiod = options.isSet("zperiod");
    const auto zmin = has_zperiod ? 0.0 : options["ZMIN"].withDefault(0.0);
    const auto zmax = has_zperiod ? 1.0 / options["zperiod"].withDefault(1.0)
                                  : options["ZMAX"].withDefault(1.0);

    const auto default_dz = (zmax - zmin) * TWOPI / nz;

    mesh->get(dz, "dz", default_dz);
  }

  // Diagonal components of metric tensor g^{ij} (default to 1)
  mesh->get(g11, "g11", 1.0);
  g11 = interpolateAndExtrapolate(g11, location, extrapolate_x, extrapolate_y);
  mesh->get(g22, "g22", 1.0);
  g22 = interpolateAndExtrapolate(g22, location, extrapolate_x, extrapolate_y);
  mesh->get(g33, "g33", 1.0);
  g33 = interpolateAndExtrapolate(g33, location, extrapolate_x, extrapolate_y);

  // Off-diagonal elements. Default to 0
  mesh->get(g12, "g12", 0.0);
  g12 = interpolateAndExtrapolate(g12, location, extrapolate_x, extrapolate_y);
  mesh->get(g13, "g13", 0.0);
  g13 = interpolateAndExtrapolate(g13, location, extrapolate_x, extrapolate_y);
  mesh->get(g23, "g23", 0.0);
  g23 = interpolateAndExtrapolate(g23, location, extrapolate_x, extrapolate_y);

  // Check input metrics
  // Diagonal metric components should be finite
  bout::checkFinite(g11, "g11", "RGN_NOCORNERS");
  bout::checkFinite(g22, "g22", "RGN_NOCORNERS");
  bout::checkFinite(g33, "g33", "RGN_NOCORNERS");
  // Diagonal metric components should be positive
  bout::checkPositive(g11, "g11", "RGN_NOCORNERS");
  bout::checkPositive(g22, "g22", "RGN_NOCORNERS");
  bout::checkPositive(g33, "g33", "RGN_NOCORNERS");
  // Off-diagonal metric components should be finite
  bout::checkFinite(g12, "g12", "RGN_NOCORNERS");
  bout::checkFinite(g13, "g13", "RGN_NOCORNERS");
  bout::checkFinite(g23, "g23", "RGN_NOCORNERS");

  /// Find covariant metric components
  auto covariant_component_names = {"g_11", "g_22", "g_33", "g_12", "g_13", "g_23"};
  auto source_has_component = [&mesh] (const std::string& name) {
    return mesh->sourceHasVar(name);
  };
  // Check if any of the components are present
  if (std::any_of(begin(covariant_component_names), end(covariant_component_names),
                  source_has_component)) {
    // Check that all components are present
    if (std::all_of(begin(covariant_component_names), end(covariant_component_names),
                    source_has_component)) {
      mesh->get(g_11, "g_11");
      mesh->get(g_22, "g_22");
      mesh->get(g_33, "g_33");
      mesh->get(g_12, "g_12");
      mesh->get(g_13, "g_13");
      mesh->get(g_23, "g_23");

      output_warn.write("\tWARNING! Covariant components of metric tensor set manually. "
                        "Contravariant components NOT recalculated\n");

    } else {
      output_warn.write("Not all covariant components of metric tensor found. "
                        "Calculating all from the contravariant tensor\n");
      /// Calculate contravariant metric components if not found
      if (calcCovariant("RGN_NOCORNERS")) {
        throw BoutException("Error in calcCovariant call");
      }
    }
  } else {
    /// Calculate contravariant metric components if not found
    if (calcCovariant("RGN_NOCORNERS")) {
      throw BoutException("Error in calcCovariant call");
    }
  }
  // More robust to extrapolate derived quantities directly, rather than
  // deriving from extrapolated covariant metric components
  g_11 = interpolateAndExtrapolate(g_11, location, extrapolate_x, extrapolate_y);
  g_22 = interpolateAndExtrapolate(g_22, location, extrapolate_x, extrapolate_y);
  g_33 = interpolateAndExtrapolate(g_33, location, extrapolate_x, extrapolate_y);
  g_12 = interpolateAndExtrapolate(g_12, location, extrapolate_x, extrapolate_y);
  g_13 = interpolateAndExtrapolate(g_13, location, extrapolate_x, extrapolate_y);
  g_23 = interpolateAndExtrapolate(g_23, location, extrapolate_x, extrapolate_y);

  // Check covariant metrics
  // Diagonal metric components should be finite
  bout::checkFinite(g_11, "g_11", "RGN_NOCORNERS");
  bout::checkFinite(g_22, "g_22", "RGN_NOCORNERS");
  bout::checkFinite(g_33, "g_33", "RGN_NOCORNERS");
  // Diagonal metric components should be positive
  bout::checkPositive(g_11, "g_11", "RGN_NOCORNERS");
  bout::checkPositive(g_22, "g_22", "RGN_NOCORNERS");
  bout::checkPositive(g_33, "g_33", "RGN_NOCORNERS");
  // Off-diagonal metric components should be finite
  bout::checkFinite(g_12, "g_12", "RGN_NOCORNERS");
  bout::checkFinite(g_13, "g_13", "RGN_NOCORNERS");
  bout::checkFinite(g_23, "g_23", "RGN_NOCORNERS");

  /// Calculate Jacobian and Bxy
  if (jacobian())
    throw BoutException("Error in jacobian call");

  // Attempt to read J from the grid file
  Field2D Jcalc = J;
  if (mesh->get(J, "J")) {
    output_warn.write(
        "\tWARNING: Jacobian 'J' not found. Calculating from metric tensor\n");
    J = Jcalc;
  } else {
    J = interpolateAndExtrapolate(J, location, extrapolate_x, extrapolate_y);

    // Compare calculated and loaded values
    output_warn.write("\tMaximum difference in J is {:e}\n", max(abs(J - Jcalc)));

    // Re-evaluate Bxy using new J
    Bxy = sqrt(g_22) / J;
  }

  // Attempt to read Bxy from the grid file
  Field2D Bcalc = Bxy;
  if (mesh->get(Bxy, "Bxy")) {
    output_warn.write("\tWARNING: Magnitude of B field 'Bxy' not found. Calculating from "
                      "metric tensor\n");
    Bxy = Bcalc;
  } else {
    Bxy = interpolateAndExtrapolate(Bxy, location, extrapolate_x, extrapolate_y);

    output_warn.write("\tMaximum difference in Bxy is {:e}\n", max(abs(Bxy - Bcalc)));
    // Check Bxy
    bout::checkFinite(Bxy, "Bxy", "RGN_NOCORNERS");
    bout::checkPositive(Bxy, "Bxy", "RGN_NOCORNERS");
  }

  //////////////////////////////////////////////////////
  /// Calculate Christoffel symbols. Needs communication
  if (geometry()) {
    throw BoutException("Differential geometry failed\n");
  }

  if (mesh->get(ShiftTorsion, "ShiftTorsion")) {
    output_warn.write(
        "\tWARNING: No Torsion specified for zShift. Derivatives may not be correct\n");
    ShiftTorsion = 0.0;
  }
  ShiftTorsion = interpolateAndExtrapolate(ShiftTorsion, location, extrapolate_x, extrapolate_y);

  //////////////////////////////////////////////////////

  if (mesh->IncIntShear) {
    if (mesh->get(IntShiftTorsion, "IntShiftTorsion")) {
      output_warn.write("\tWARNING: No Integrated torsion specified\n");
      IntShiftTorsion = 0.0;
    }
    IntShiftTorsion = interpolateAndExtrapolate(IntShiftTorsion, location, extrapolate_x, extrapolate_y);
  } else {
    // IntShiftTorsion will not be used, but set to zero to avoid uninitialized field
    IntShiftTorsion = 0.;
  }

  setParallelTransform(options);
}

Coordinates::Coordinates(Mesh* mesh, Options* options, const CELL_LOC loc,
      const Coordinates* coords_in, bool force_interpolate_from_centre)
    : dx(1, mesh), dy(1, mesh), dz(1), d1_dx(mesh), d1_dy(mesh), J(1, mesh), Bxy(1, mesh),
      // Identity metric tensor
      g11(1, mesh), g22(1, mesh), g33(1, mesh), g12(0, mesh), g13(0, mesh), g23(0, mesh),
      g_11(1, mesh), g_22(1, mesh), g_33(1, mesh), g_12(0, mesh), g_13(0, mesh),
      g_23(0, mesh), G1_11(mesh), G1_22(mesh), G1_33(mesh), G1_12(mesh), G1_13(mesh),
      G1_23(mesh), G2_11(mesh), G2_22(mesh), G2_33(mesh), G2_12(mesh), G2_13(mesh),
      G2_23(mesh), G3_11(mesh), G3_22(mesh), G3_33(mesh), G3_12(mesh), G3_13(mesh),
      G3_23(mesh), G1(mesh), G2(mesh), G3(mesh), ShiftTorsion(mesh),
      IntShiftTorsion(mesh), localmesh(mesh), location(loc) {

  std::string suffix = getLocationSuffix(location);

  nz = mesh->LocalNz;

  dz = coords_in->dz;

  // Default to true in case staggered quantities are not read from file
  bool extrapolate_x = true;
  bool extrapolate_y = true;

  if (!force_interpolate_from_centre && mesh->sourceHasVar("dx"+suffix)) {

    extrapolate_x = not mesh->sourceHasXBoundaryGuards();
    extrapolate_y = not mesh->sourceHasYBoundaryGuards();

    if (extrapolate_x) {
      output_warn.write(_("WARNING: extrapolating input mesh quantities into x-boundary "
            "cells\n"));
    }

    if (extrapolate_y) {
      output_warn.write(_("WARNING: extrapolating input mesh quantities into y-boundary "
            "cells\n"));
    }

    getAtLoc(mesh, dx, "dx", suffix, location, 1.0);
    dx = interpolateAndExtrapolate(dx, location, extrapolate_x, extrapolate_y);

    if (mesh->periodicX) {
      mesh->communicate(dx);
    }

    getAtLoc(mesh, dy, "dy", suffix, location, 1.0);
    dy = interpolateAndExtrapolate(dy, location, extrapolate_x, extrapolate_y);

    // grid data source has staggered fields, so read instead of interpolating
    // Diagonal components of metric tensor g^{ij} (default to 1)
    getAtLoc(mesh, g11, "g11", suffix, location, 1.0);
    g11 = interpolateAndExtrapolate(g11, location, extrapolate_x, extrapolate_y);
    getAtLoc(mesh, g22, "g22", suffix, location, 1.0);
    g22 = interpolateAndExtrapolate(g22, location, extrapolate_x, extrapolate_y);
    getAtLoc(mesh, g33, "g33", suffix, location, 1.0);
    g33 = interpolateAndExtrapolate(g33, location, extrapolate_x, extrapolate_y);
    getAtLoc(mesh, g12, "g12", suffix, location, 0.0);
    g12 = interpolateAndExtrapolate(g12, location, extrapolate_x, extrapolate_y);
    getAtLoc(mesh, g13, "g13", suffix, location, 0.0);
    g13 = interpolateAndExtrapolate(g13, location, extrapolate_x, extrapolate_y);
    getAtLoc(mesh, g23, "g23", suffix, location, 0.0);
    g23 = interpolateAndExtrapolate(g23, location, extrapolate_x, extrapolate_y);

    /// Find covariant metric components
    auto covariant_component_names = {"g_11", "g_22", "g_33", "g_12", "g_13", "g_23"};
    auto source_has_component = [&suffix, &mesh] (const std::string& name) {
      return mesh->sourceHasVar(name + suffix);
    };
    // Check if any of the components are present
    if (std::any_of(begin(covariant_component_names), end(covariant_component_names),
                    source_has_component)) {
      // Check that all components are present
      if (std::all_of(begin(covariant_component_names), end(covariant_component_names),
                      source_has_component)) {

        getAtLoc(mesh, g_11, "g_11", suffix, location);
        getAtLoc(mesh, g_22, "g_22", suffix, location);
        getAtLoc(mesh, g_33, "g_33", suffix, location);
        getAtLoc(mesh, g_12, "g_12", suffix, location);
        getAtLoc(mesh, g_13, "g_13", suffix, location);
        getAtLoc(mesh, g_23, "g_23", suffix, location);

        output_warn.write("\tWARNING! Staggered covariant components of metric tensor set manually. "
                          "Contravariant components NOT recalculated\n");

      } else {
        output_warn.write("Not all staggered covariant components of metric tensor found. "
                          "Calculating all from the contravariant tensor\n");
        /// Calculate contravariant metric components if not found
        if (calcCovariant("RGN_NOCORNERS")) {
          throw BoutException("Error in staggered calcCovariant call");
        }
      }
    } else {
      /// Calculate contravariant metric components if not found
      if (calcCovariant("RGN_NOCORNERS")) {
        throw BoutException("Error in staggered calcCovariant call");
      }
    }
    // More robust to extrapolate derived quantities directly, rather than
    // deriving from extrapolated covariant metric components
    g_11 = interpolateAndExtrapolate(g_11, location, extrapolate_x, extrapolate_y);
    g_22 = interpolateAndExtrapolate(g_22, location, extrapolate_x, extrapolate_y);
    g_33 = interpolateAndExtrapolate(g_33, location, extrapolate_x, extrapolate_y);
    g_12 = interpolateAndExtrapolate(g_12, location, extrapolate_x, extrapolate_y);
    g_13 = interpolateAndExtrapolate(g_13, location, extrapolate_x, extrapolate_y);
    g_23 = interpolateAndExtrapolate(g_23, location, extrapolate_x, extrapolate_y);

    /// Calculate Jacobian and Bxy
    if (jacobian()) {
      throw BoutException("Error in jacobian call while constructing staggered Coordinates");
    }

    checkStaggeredGet(mesh, "ShiftTorsion", suffix);
    if (mesh->get(ShiftTorsion, "ShiftTorsion"+suffix)) {
      output_warn.write("\tWARNING: No Torsion specified for zShift. Derivatives may not be correct\n");
      ShiftTorsion = 0.0;
    }
    ShiftTorsion.setLocation(location);
    ShiftTorsion = interpolateAndExtrapolate(ShiftTorsion, location, extrapolate_x, extrapolate_y);

    //////////////////////////////////////////////////////

    if (mesh->IncIntShear) {
      checkStaggeredGet(mesh, "IntShiftTorsion", suffix);
      if (mesh->get(IntShiftTorsion, "IntShiftTorsion"+suffix)) {
        output_warn.write("\tWARNING: No Integrated torsion specified\n");
        IntShiftTorsion = 0.0;
      }
      IntShiftTorsion.setLocation(location);
      IntShiftTorsion = interpolateAndExtrapolate(IntShiftTorsion, location, extrapolate_x, extrapolate_y);
    } else {
      // IntShiftTorsion will not be used, but set to zero to avoid uninitialized field
      IntShiftTorsion = 0.;
    }
  } else {
    // Interpolate fields from coords_in

    dx = interpolateAndExtrapolate(coords_in->dx, location);
    dy = interpolateAndExtrapolate(coords_in->dy, location);

    // Diagonal components of metric tensor g^{ij}
    g11 = interpolateAndExtrapolate(coords_in->g11, location);
    g22 = interpolateAndExtrapolate(coords_in->g22, location);
    g33 = interpolateAndExtrapolate(coords_in->g33, location);

    // Off-diagonal elements.
    g12 = interpolateAndExtrapolate(coords_in->g12, location);
    g13 = interpolateAndExtrapolate(coords_in->g13, location);
    g23 = interpolateAndExtrapolate(coords_in->g23, location);

    // 3x3 matrix inversion can exaggerate small interpolation errors, so it is
    // more robust to interpolate and extrapolate derived quantities directly,
    // rather than deriving from interpolated/extrapolated covariant metric
    // components
    g_11 = interpolateAndExtrapolate(coords_in->g_11, location);
    g_22 = interpolateAndExtrapolate(coords_in->g_22, location);
    g_33 = interpolateAndExtrapolate(coords_in->g_33, location);
    g_12 = interpolateAndExtrapolate(coords_in->g_12, location);
    g_13 = interpolateAndExtrapolate(coords_in->g_13, location);
    g_23 = interpolateAndExtrapolate(coords_in->g_23, location);

    J = interpolateAndExtrapolate(coords_in->J, location);
    Bxy = interpolateAndExtrapolate(coords_in->Bxy, location);

    bout::checkFinite(J, "The Jacobian", "RGN_NOCORNERS");
    bout::checkPositive(J, "The Jacobian", "RGN_NOCORNERS");
    bout::checkFinite(Bxy, "Bxy", "RGN_NOCORNERS");
    bout::checkPositive(Bxy, "Bxy", "RGN_NOCORNERS");

    ShiftTorsion = interpolateAndExtrapolate(coords_in->ShiftTorsion, location);

    if (mesh->IncIntShear) {
      IntShiftTorsion = interpolateAndExtrapolate(coords_in->IntShiftTorsion, location);
    }
  }

  // Check input metrics
  // Diagonal metric components should be finite
  bout::checkFinite(g11, "g11", "RGN_NOCORNERS");
  bout::checkFinite(g22, "g22", "RGN_NOCORNERS");
  bout::checkFinite(g33, "g33", "RGN_NOCORNERS");
  bout::checkFinite(g_11, "g_11", "RGN_NOCORNERS");
  bout::checkFinite(g_22, "g_22", "RGN_NOCORNERS");
  bout::checkFinite(g_33, "g_33", "RGN_NOCORNERS");
  // Diagonal metric components should be positive
  bout::checkPositive(g11, "g11", "RGN_NOCORNERS");
  bout::checkPositive(g22, "g22", "RGN_NOCORNERS");
  bout::checkPositive(g33, "g33", "RGN_NOCORNERS");
  bout::checkPositive(g_11, "g_11", "RGN_NOCORNERS");
  bout::checkPositive(g_22, "g_22", "RGN_NOCORNERS");
  bout::checkPositive(g_33, "g_33", "RGN_NOCORNERS");
  // Off-diagonal metric components should be finite
  bout::checkFinite(g12, "g12", "RGN_NOCORNERS");
  bout::checkFinite(g13, "g13", "RGN_NOCORNERS");
  bout::checkFinite(g23, "g23", "RGN_NOCORNERS");
  bout::checkFinite(g_12, "g_12", "RGN_NOCORNERS");
  bout::checkFinite(g_13, "g_13", "RGN_NOCORNERS");
  bout::checkFinite(g_23, "g_23", "RGN_NOCORNERS");

  //////////////////////////////////////////////////////
  /// Calculate Christoffel symbols. Needs communication
  if (geometry(false, force_interpolate_from_centre)) {
    throw BoutException("Differential geometry failed while constructing staggered Coordinates");
  }

  setParallelTransform(options);
}

void Coordinates::outputVars(Datafile& file) {
  const std::string loc_string = (location == CELL_CENTRE) ? "" : "_"+toString(location);

  file.addOnce(dx, "dx" + loc_string);
  file.addOnce(dy, "dy" + loc_string);
  file.addOnce(dz, "dz" + loc_string);

  file.addOnce(g11, "g11" + loc_string);
  file.addOnce(g22, "g22" + loc_string);
  file.addOnce(g33, "g33" + loc_string);
  file.addOnce(g12, "g12" + loc_string);
  file.addOnce(g13, "g13" + loc_string);
  file.addOnce(g23, "g23" + loc_string);

  file.addOnce(g_11, "g_11" + loc_string);
  file.addOnce(g_22, "g_22" + loc_string);
  file.addOnce(g_33, "g_33" + loc_string);
  file.addOnce(g_12, "g_12" + loc_string);
  file.addOnce(g_13, "g_13" + loc_string);
  file.addOnce(g_23, "g_23" + loc_string);

  file.addOnce(J, "J" + loc_string);
  file.addOnce(Bxy, "Bxy" + loc_string);

  file.addOnce(G1, "G1" + loc_string);
  file.addOnce(G2, "G2" + loc_string);
  file.addOnce(G3, "G3" + loc_string);

  getParallelTransform().outputVars(file);
}

int Coordinates::geometry(bool recalculate_staggered,
    bool force_interpolate_from_centre) {
  TRACE("Coordinates::geometry");

  output_progress.write("Calculating differential geometry terms\n");

  if (min(abs(dx)) < 1e-8)
    throw BoutException("dx magnitude less than 1e-8");

  if (min(abs(dy)) < 1e-8)
    throw BoutException("dy magnitude less than 1e-8");

  if (fabs(dz) < 1e-8)
    throw BoutException("dz magnitude less than 1e-8");

  // Check input metrics
  // Diagonal metric components should be finite
  bout::checkFinite(g11, "g11", "RGN_NOCORNERS");
  bout::checkFinite(g22, "g22", "RGN_NOCORNERS");
  bout::checkFinite(g33, "g33", "RGN_NOCORNERS");
  // Diagonal metric components should be positive
  bout::checkPositive(g11, "g11", "RGN_NOCORNERS");
  bout::checkPositive(g22, "g22", "RGN_NOCORNERS");
  bout::checkPositive(g33, "g33", "RGN_NOCORNERS");
  // Off-diagonal metric components should be finite
  bout::checkFinite(g12, "g12", "RGN_NOCORNERS");
  bout::checkFinite(g13, "g13", "RGN_NOCORNERS");
  bout::checkFinite(g23, "g23", "RGN_NOCORNERS");

  // Diagonal metric components should be finite
  bout::checkFinite(g_11, "g_11", "RGN_NOCORNERS");
  bout::checkFinite(g_22, "g_22", "RGN_NOCORNERS");
  bout::checkFinite(g_33, "g_33", "RGN_NOCORNERS");
  // Diagonal metric components should be positive
  bout::checkPositive(g_11, "g_11", "RGN_NOCORNERS");
  bout::checkPositive(g_22, "g_22", "RGN_NOCORNERS");
  bout::checkPositive(g_33, "g_33", "RGN_NOCORNERS");
  // Off-diagonal metric components should be finite
  bout::checkFinite(g_12, "g_12", "RGN_NOCORNERS");
  bout::checkFinite(g_13, "g_13", "RGN_NOCORNERS");
  bout::checkFinite(g_23, "g_23", "RGN_NOCORNERS");

  // Calculate Christoffel symbol terms (18 independent values)
  // Note: This calculation is completely general: metric
  // tensor can be 2D or 3D. For 2D, all DDZ terms are zero

  G1_11 = 0.5 * g11 * DDX(g_11) + g12 * (DDX(g_12) - 0.5 * DDY(g_11))
          + g13 * (DDX(g_13) - 0.5 * DDZ(g_11));
  G1_22 = g11 * (DDY(g_12) - 0.5 * DDX(g_22)) + 0.5 * g12 * DDY(g_22)
          + g13 * (DDY(g_23) - 0.5 * DDZ(g_22));
  G1_33 = g11 * (DDZ(g_13) - 0.5 * DDX(g_33)) + g12 * (DDZ(g_23) - 0.5 * DDY(g_33))
          + 0.5 * g13 * DDZ(g_33);
  G1_12 = 0.5 * g11 * DDY(g_11) + 0.5 * g12 * DDX(g_22)
          + 0.5 * g13 * (DDY(g_13) + DDX(g_23) - DDZ(g_12));
  G1_13 = 0.5 * g11 * DDZ(g_11) + 0.5 * g12 * (DDZ(g_12) + DDX(g_23) - DDY(g_13))
          + 0.5 * g13 * DDX(g_33);
  G1_23 = 0.5 * g11 * (DDZ(g_12) + DDY(g_13) - DDX(g_23))
          + 0.5 * g12 * (DDZ(g_22) + DDY(g_23) - DDY(g_23))
          // + 0.5 *g13*(DDZ(g_32) + DDY(g_33) - DDZ(g_23));
          // which equals
          + 0.5 * g13 * DDY(g_33);

  G2_11 = 0.5 * g12 * DDX(g_11) + g22 * (DDX(g_12) - 0.5 * DDY(g_11))
          + g23 * (DDX(g_13) - 0.5 * DDZ(g_11));
  G2_22 = g12 * (DDY(g_12) - 0.5 * DDX(g_22)) + 0.5 * g22 * DDY(g_22)
          + g23 * (DDY(g23) - 0.5 * DDZ(g_22));
  G2_33 = g12 * (DDZ(g_13) - 0.5 * DDX(g_33)) + g22 * (DDZ(g_23) - 0.5 * DDY(g_33))
          + 0.5 * g23 * DDZ(g_33);
  G2_12 = 0.5 * g12 * DDY(g_11) + 0.5 * g22 * DDX(g_22)
          + 0.5 * g23 * (DDY(g_13) + DDX(g_23) - DDZ(g_12));
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
  G2_23 = 0.5 * g12 * (DDZ(g_12) + DDY(g_13) - DDX(g_23)) + 0.5 * g22 * DDZ(g_22)
          + 0.5 * g23 * DDY(g_33);

  G3_11 = 0.5 * g13 * DDX(g_11) + g23 * (DDX(g_12) - 0.5 * DDY(g_11))
          + g33 * (DDX(g_13) - 0.5 * DDZ(g_11));
  G3_22 = g13 * (DDY(g_12) - 0.5 * DDX(g_22)) + 0.5 * g23 * DDY(g_22)
          + g33 * (DDY(g_23) - 0.5 * DDZ(g_22));
  G3_33 = g13 * (DDZ(g_13) - 0.5 * DDX(g_33)) + g23 * (DDZ(g_23) - 0.5 * DDY(g_33))
          + 0.5 * g33 * DDZ(g_33);
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
  G3_13 = 0.5 * g13 * DDZ(g_11) + 0.5 * g23 * (DDZ(g_12) + DDX(g_23) - DDY(g_13))
          + 0.5 * g33 * DDX(g_33);
  G3_23 = 0.5 * g13 * (DDZ(g_12) + DDY(g_13) - DDX(g_23)) + 0.5 * g23 * DDZ(g_22)
          + 0.5 * g33 * DDY(g_33);

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

  // Set boundary guard cells of Christoffel symbol terms
  // Ideally, when location is staggered, we would set the upper/outer boundary point
  // correctly rather than by extrapolating here: e.g. if location==CELL_YLOW and we are
  // at the upper y-boundary the x- and z-derivatives at yend+1 at the boundary can be
  // calculated because the guard cells are available, while the y-derivative could be
  // calculated from the CELL_CENTRE metric components (which have guard cells available
  // past the boundary location). This would avoid the problem that the y-boundary on the
  // CELL_YLOW grid is at a 'guard cell' location (yend+1).
  // However, the above would require lots of special handling, so just extrapolate for
  // now.
  G1_11 = interpolateAndExtrapolate(G1_11, location, true, true, true);
  G1_22 = interpolateAndExtrapolate(G1_22, location, true, true, true);
  G1_33 = interpolateAndExtrapolate(G1_33, location, true, true, true);
  G1_12 = interpolateAndExtrapolate(G1_12, location, true, true, true);
  G1_13 = interpolateAndExtrapolate(G1_13, location, true, true, true);
  G1_23 = interpolateAndExtrapolate(G1_23, location, true, true, true);

  G2_11 = interpolateAndExtrapolate(G2_11, location, true, true, true);
  G2_22 = interpolateAndExtrapolate(G2_22, location, true, true, true);
  G2_33 = interpolateAndExtrapolate(G2_33, location, true, true, true);
  G2_12 = interpolateAndExtrapolate(G2_12, location, true, true, true);
  G2_13 = interpolateAndExtrapolate(G2_13, location, true, true, true);
  G2_23 = interpolateAndExtrapolate(G2_23, location, true, true, true);

  G3_11 = interpolateAndExtrapolate(G3_11, location, true, true, true);
  G3_22 = interpolateAndExtrapolate(G3_22, location, true, true, true);
  G3_33 = interpolateAndExtrapolate(G3_33, location, true, true, true);
  G3_12 = interpolateAndExtrapolate(G3_12, location, true, true, true);
  G3_13 = interpolateAndExtrapolate(G3_13, location, true, true, true);
  G3_23 = interpolateAndExtrapolate(G3_23, location, true, true, true);

  G1 = interpolateAndExtrapolate(G1, location, true, true, true);
  G2 = interpolateAndExtrapolate(G2, location, true, true, true);
  G3 = interpolateAndExtrapolate(G3, location, true, true, true);

  //////////////////////////////////////////////////////
  /// Non-uniform meshes. Need to use DDX, DDY

  OPTION(Options::getRoot(), non_uniform, true);

  Field2D d2x(localmesh), d2y(localmesh); // d^2 x / d i^2
  // Read correction for non-uniform meshes
  std::string suffix = getLocationSuffix(location);
  if (location == CELL_CENTRE or (!force_interpolate_from_centre
                      and localmesh->sourceHasVar("dx"+suffix))) {
    bool extrapolate_x = not localmesh->sourceHasXBoundaryGuards();
    bool extrapolate_y = not localmesh->sourceHasYBoundaryGuards();

    if (localmesh->get(d2x, "d2x"+suffix)) {
      output_warn.write(
          "\tWARNING: differencing quantity 'd2x' not found. Calculating from dx\n");
      d1_dx = bout::derivatives::index::DDX(1. / dx); // d/di(1/dx)

      localmesh->communicate(d1_dx);
      d1_dx = interpolateAndExtrapolate(d1_dx, location, true, true, true);
    } else {
      d2x.setLocation(location);
      // set boundary cells if necessary
      d2x = interpolateAndExtrapolate(d2x, location, extrapolate_x, extrapolate_y);

      d1_dx = -d2x / (dx * dx);
    }

    if (localmesh->get(d2y, "d2y"+suffix)) {
      output_warn.write(
          "\tWARNING: differencing quantity 'd2y' not found. Calculating from dy\n");
      d1_dy = bout::derivatives::index::DDY(1. / dy); // d/di(1/dy)

      localmesh->communicate(d1_dy);
      d1_dy = interpolateAndExtrapolate(d1_dy, location, true, true, true);
    } else {
      d2y.setLocation(location);
      // set boundary cells if necessary
      d2y = interpolateAndExtrapolate(d2y, location, extrapolate_x, extrapolate_y);

      d1_dy = -d2y / (dy * dy);
    }
  } else {
    if (localmesh->get(d2x, "d2x")) {
      output_warn.write(
          "\tWARNING: differencing quantity 'd2x' not found. Calculating from dx\n");
      d1_dx = bout::derivatives::index::DDX(1. / dx); // d/di(1/dx)

      localmesh->communicate(d1_dx);
      d1_dx = interpolateAndExtrapolate(d1_dx, location, true, true, true);
    } else {
      // Shift d2x to our location
      d2x = interpolateAndExtrapolate(d2x, location);

      d1_dx = -d2x / (dx * dx);
    }

    if (localmesh->get(d2y, "d2y")) {
      output_warn.write(
          "\tWARNING: differencing quantity 'd2y' not found. Calculating from dy\n");
      d1_dy = bout::derivatives::index::DDY(1. / dy); // d/di(1/dy)

      localmesh->communicate(d1_dy);
      d1_dy = interpolateAndExtrapolate(d1_dy, location, true, true, true);
    } else {
      // Shift d2y to our location
      d2y = interpolateAndExtrapolate(d2y, location);

      d1_dy = -d2y / (dy * dy);
    }
  }

  if (location == CELL_CENTRE && recalculate_staggered) {
    // Re-calculate interpolated Coordinates at staggered locations
    localmesh->recalculateStaggeredCoordinates();
  }

  return 0;
}

int Coordinates::calcCovariant(const std::string& region) {
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

  BOUT_FOR_SERIAL(i, g11.getRegion(region)) {
    a(0, 0) = g11[i];
    a(1, 1) = g22[i];
    a(2, 2) = g33[i];

    a(0, 1) = a(1, 0) = g12[i];
    a(1, 2) = a(2, 1) = g23[i];
    a(0, 2) = a(2, 0) = g13[i];

    if (invert3x3(a)) {
      output_error.write("\tERROR: metric tensor is singular at ({:d}, {:d})\n", i.x(), i.y());
      return 1;
    }

    g_11[i] = a(0, 0);
    g_22[i] = a(1, 1);
    g_33[i] = a(2, 2);

    g_12[i] = a(0, 1);
    g_13[i] = a(0, 2);
    g_23[i] = a(1, 2);
  }

  BoutReal maxerr;
  maxerr = BOUTMAX(max(abs((g_11 * g11 + g_12 * g12 + g_13 * g13) - 1)),
                   max(abs((g_12 * g12 + g_22 * g22 + g_23 * g23) - 1)),
                   max(abs((g_13 * g13 + g_23 * g23 + g_33 * g33) - 1)));

  output_info.write("\tLocal maximum error in diagonal inversion is {:e}\n", maxerr);

  maxerr = BOUTMAX(max(abs(g_11 * g12 + g_12 * g22 + g_13 * g23)),
                   max(abs(g_11 * g13 + g_12 * g23 + g_13 * g33)),
                   max(abs(g_12 * g13 + g_22 * g23 + g_23 * g33)));

  output_info.write("\tLocal maximum error in off-diagonal inversion is {:e}\n", maxerr);

  return 0;
}

int Coordinates::calcContravariant(const std::string& region) {
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

  BOUT_FOR_SERIAL(i, g_11.getRegion(region)) {
    a(0, 0) = g_11[i];
    a(1, 1) = g_22[i];
    a(2, 2) = g_33[i];

    a(0, 1) = a(1, 0) = g_12[i];
    a(1, 2) = a(2, 1) = g_23[i];
    a(0, 2) = a(2, 0) = g_13[i];

    if (invert3x3(a)) {
      output_error.write("\tERROR: metric tensor is singular at ({:d}, {:d})\n", i.x(), i.y());
      return 1;
    }

    g11[i] = a(0, 0);
    g22[i] = a(1, 1);
    g33[i] = a(2, 2);

    g12[i] = a(0, 1);
    g13[i] = a(0, 2);
    g23[i] = a(1, 2);
  }

  BoutReal maxerr;
  maxerr = BOUTMAX(max(abs((g_11 * g11 + g_12 * g12 + g_13 * g13) - 1)),
                   max(abs((g_12 * g12 + g_22 * g22 + g_23 * g23) - 1)),
                   max(abs((g_13 * g13 + g_23 * g23 + g_33 * g33) - 1)));

  output_info.write("\tMaximum error in diagonal inversion is {:e}\n", maxerr);

  maxerr = BOUTMAX(max(abs(g_11 * g12 + g_12 * g22 + g_13 * g23)),
                   max(abs(g_11 * g13 + g_12 * g23 + g_13 * g33)),
                   max(abs(g_12 * g13 + g_22 * g23 + g_23 * g33)));

  output_info.write("\tMaximum error in off-diagonal inversion is {:e}\n", maxerr);
  return 0;
}

int Coordinates::jacobian() {
  TRACE("Coordinates::jacobian");
  // calculate Jacobian using g^-1 = det[g^ij], J = sqrt(g)

  const bool extrapolate_x = not localmesh->sourceHasXBoundaryGuards();
  const bool extrapolate_y = not localmesh->sourceHasYBoundaryGuards();

  Field2D g = g11 * g22 * g33 + 2.0 * g12 * g13 * g23 - g11 * g23 * g23 - g22 * g13 * g13
              - g33 * g12 * g12;

  // Check that g is positive
  bout::checkPositive(g, "The determinant of g^ij", "RGN_NOBNDRY");

  J = 1. / sqrt(g);
  // More robust to extrapolate derived quantities directly, rather than
  // deriving from extrapolated covariant metric components
  J = interpolateAndExtrapolate(J, location, extrapolate_x, extrapolate_y);

  // Check jacobian
  bout::checkFinite(J, "The Jacobian", "RGN_NOCORNERS");
  bout::checkPositive(J, "The Jacobian", "RGN_NOCORNERS");
  if (min(abs(J)) < 1.0e-10) {
    throw BoutException("\tERROR: Jacobian becomes very small\n");
  }

  bout::checkPositive(g_22, "g_22", "RGN_NOCORNERS");

  Bxy = sqrt(g_22) / J;
  Bxy = interpolateAndExtrapolate(Bxy, location, extrapolate_x, extrapolate_y);

  bout::checkFinite(Bxy, "Bxy", "RGN_NOCORNERS");
  bout::checkPositive(Bxy, "Bxy", "RGN_NOCORNERS");

  return 0;
}

namespace {
// Utility function for fixing up guard cells of zShift
void fixZShiftGuards(Field2D& zShift) {
  auto localmesh = zShift.getMesh();

  // extrapolate into boundary guard cells if necessary
  zShift = interpolateAndExtrapolate(zShift, zShift.getLocation(),
      not localmesh->sourceHasXBoundaryGuards(),
      not localmesh->sourceHasYBoundaryGuards());

  // make sure zShift has been communicated
  localmesh->communicate(zShift);

  // Correct guard cells for discontinuity of zShift at poloidal branch cut
  for (int x = 0; x < localmesh->LocalNx; x++) {
    const auto lower = localmesh->hasBranchCutLower(x);
    if (lower.first) {
      for (int y = 0; y < localmesh->ystart; y++) {
        zShift(x, y) -= lower.second;
      }
    }
    const auto upper = localmesh->hasBranchCutUpper(x);
    if (upper.first) {
      for (int y = localmesh->yend + 1; y < localmesh->LocalNy; y++) {
        zShift(x, y) += upper.second;
      }
    }
  }
}
}

void Coordinates::setParallelTransform(Options* options) {

  auto ptoptions = options->getSection("paralleltransform");

  std::string ptstr;
  ptoptions->get("type", ptstr, "identity");

  // Convert to lower case for comparison
  ptstr = lowercase(ptstr);

  if(ptstr == "identity") {
    // Identity method i.e. no transform needed
    transform = bout::utils::make_unique<ParallelTransformIdentity>(*localmesh,
                                                                    ptoptions);

  } else if (ptstr == "shifted") {
    // Shifted metric method

    Field2D zShift{localmesh};

    // Read the zShift angle from the mesh
    std::string suffix = getLocationSuffix(location);
    if (localmesh->sourceHasVar("dx"+suffix)) {
      // Grid file has variables at this location, so should be able to read
      checkStaggeredGet(localmesh, "zShift", suffix);
      if (localmesh->get(zShift, "zShift"+suffix)) {
        // No zShift variable. Try qinty in BOUT grid files
        if (localmesh->get(zShift, "qinty"+suffix)) {
          // Failed to find either variable, cannot use ShiftedMetric
          throw BoutException("Could not read zShift"+suffix+" from grid file");
        }
      }
      zShift.setLocation(location);
    } else {
      Field2D zShift_centre;
      if (localmesh->get(zShift_centre, "zShift")) {
        // No zShift variable. Try qinty in BOUT grid files
        if (localmesh->get(zShift_centre, "qinty")) {
          // Failed to find either variable, cannot use ShiftedMetric
          throw BoutException("Could not read zShift"+suffix+" from grid file");
        }
      }

      fixZShiftGuards(zShift_centre);

      zShift = interpolateAndExtrapolate(zShift_centre, location);
    }

    fixZShiftGuards(zShift);

    transform = bout::utils::make_unique<ShiftedMetric>(*localmesh, location, zShift,
                                                        zlength(), ptoptions);

  } else if (ptstr == "fci") {

    if (location != CELL_CENTRE) {
      throw BoutException("FCITransform is not available on staggered grids.");
    }

    // Flux Coordinate Independent method
    const bool fci_zperiodic = (*ptoptions)["z_periodic"].withDefault(true);
    transform = bout::utils::make_unique<FCITransform>(*localmesh, fci_zperiodic,
                                                       ptoptions);

  } else {
    throw BoutException(_("Unrecognised paralleltransform option.\n"
                          "Valid choices are 'identity', 'shifted', 'fci'"));
  }
}

/*******************************************************************************
 * Operators
 *
 *******************************************************************************/

Field2D Coordinates::DDX(const Field2D& f, CELL_LOC loc, const std::string& method,
    const std::string& region) {
  ASSERT1(location == loc || loc == CELL_DEFAULT);
  return bout::derivatives::index::DDX(f, loc, method, region) / dx;
}

Field2D Coordinates::DDY(const Field2D& f, CELL_LOC loc, const std::string& method,
    const std::string& region) {
  ASSERT1(location == loc || loc == CELL_DEFAULT);
  return bout::derivatives::index::DDY(f, loc, method, region) / dy;
}

Field2D Coordinates::DDZ(MAYBE_UNUSED(const Field2D& f), CELL_LOC loc,
    const std::string& UNUSED(method), const std::string& UNUSED(region)) {
  ASSERT1(location == loc || loc == CELL_DEFAULT);
  ASSERT1(f.getMesh() == localmesh);
  if (loc == CELL_DEFAULT) {
    loc = f.getLocation();
  }
  return zeroFrom(f).setLocation(loc);
}

#include <derivs.hxx>

/////////////////////////////////////////////////////////
// Parallel gradient

Field2D Coordinates::Grad_par(const Field2D& var, MAYBE_UNUSED(CELL_LOC outloc),
    const std::string& UNUSED(method)) {
  TRACE("Coordinates::Grad_par( Field2D )");
  ASSERT1(location == outloc
          || (outloc == CELL_DEFAULT && location == var.getLocation()));

  return DDY(var) / sqrt(g_22);
}

Field3D Coordinates::Grad_par(const Field3D& var, CELL_LOC outloc,
    const std::string& method) {
  TRACE("Coordinates::Grad_par( Field3D )");
  ASSERT1(location == outloc || outloc == CELL_DEFAULT);

  return ::DDY(var, outloc, method) / sqrt(g_22);
}

/////////////////////////////////////////////////////////
// Vpar_Grad_par
// vparallel times the parallel derivative along unperturbed B-field

Field2D Coordinates::Vpar_Grad_par(const Field2D& v, const Field2D& f,
    MAYBE_UNUSED(CELL_LOC outloc), const std::string& UNUSED(method)) {
  ASSERT1(location == outloc || (outloc == CELL_DEFAULT && location == f.getLocation()));
  return VDDY(v, f) / sqrt(g_22);
}

Field3D Coordinates::Vpar_Grad_par(const Field3D& v, const Field3D& f, CELL_LOC outloc,
    const std::string& method) {
  ASSERT1(location == outloc || outloc == CELL_DEFAULT);
  return VDDY(v, f, outloc, method) / sqrt(g_22);
}

/////////////////////////////////////////////////////////
// Parallel divergence

Field2D Coordinates::Div_par(const Field2D& f, CELL_LOC outloc,
    const std::string& method) {
  TRACE("Coordinates::Div_par( Field2D )");
  ASSERT1(location == outloc || outloc == CELL_DEFAULT);

  // Need Bxy at location of f, which might be different from location of this
  // Coordinates object
  Field2D Bxy_floc = f.getCoordinates()->Bxy;

  return Bxy * Grad_par(f / Bxy_floc, outloc, method);
}

Field3D Coordinates::Div_par(const Field3D& f, CELL_LOC outloc,
    const std::string& method) {
  TRACE("Coordinates::Div_par( Field3D )");
  ASSERT1(location == outloc || outloc == CELL_DEFAULT);

  // Need Bxy at location of f, which might be different from location of this
  // Coordinates object
  Field2D Bxy_floc = f.getCoordinates()->Bxy;

  if (!f.hasParallelSlices()) {
    // No yup/ydown fields. The Grad_par operator will
    // shift to field aligned coordinates
    return Bxy * Grad_par(f / Bxy_floc, outloc, method);
  }

  // Need to modify yup and ydown fields
  Field3D f_B = f / Bxy_floc;
  f_B.splitParallelSlices();
  f_B.yup() = f.yup() / Bxy_floc;
  f_B.ydown() = f.ydown() / Bxy_floc;
  return Bxy * Grad_par(f_B, outloc, method);
}

/////////////////////////////////////////////////////////
// second parallel derivative (b dot Grad)(b dot Grad)
// Note: For parallel Laplacian use Laplace_par

Field2D Coordinates::Grad2_par2(const Field2D& f, CELL_LOC outloc,
    const std::string& method) {
  TRACE("Coordinates::Grad2_par2( Field2D )");
  ASSERT1(location == outloc || (outloc == CELL_DEFAULT && location == f.getLocation()));

  Field2D sg = sqrt(g_22);
  Field2D result = DDY(1. / sg, outloc, method) * DDY(f, outloc, method) / sg
                   + D2DY2(f, outloc, method) / g_22;

  return result;
}

Field3D Coordinates::Grad2_par2(const Field3D& f, CELL_LOC outloc,
    const std::string& method) {
  TRACE("Coordinates::Grad2_par2( Field3D )");
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  ASSERT1(location == outloc);

  Field2D sg = sqrt(g_22);
  sg = DDY(1. / sg, outloc, method) / sg;

  Field3D result = ::DDY(f, outloc, method);

  Field3D r2 = D2DY2(f, outloc, method) / g_22;

  result = sg * result + r2;

  ASSERT2(result.getLocation() == outloc);

  return result;
}

/////////////////////////////////////////////////////////
// perpendicular Laplacian operator

#include <invert_laplace.hxx> // Delp2 uses same coefficients as inversion code

Field2D Coordinates::Delp2(const Field2D& f, CELL_LOC outloc, bool UNUSED(useFFT)) {
  TRACE("Coordinates::Delp2( Field2D )");
  ASSERT1(location == outloc || outloc == CELL_DEFAULT);

  Field2D result = G1 * DDX(f, outloc) + g11 * D2DX2(f, outloc);

  return result;
}

Field3D Coordinates::Delp2(const Field3D& f, CELL_LOC outloc, bool useFFT) {
  TRACE("Coordinates::Delp2( Field3D )");

  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }

  ASSERT1(location == outloc);
  ASSERT1(f.getLocation() == outloc);

  if (localmesh->GlobalNx == 1 && localmesh->GlobalNz == 1) {
    // copy mesh, location, etc
    return f * 0;
  }
  ASSERT2(localmesh->xstart > 0); // Need at least one guard cell

  Field3D result{emptyFrom(f).setLocation(outloc)};

  if (useFFT) {
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

        // No smoothing in the x direction
        for (int jx = localmesh->xstart; jx <= localmesh->xend; jx++) {
          // Perform x derivative

          dcomplex a, b, c;
          laplace_tridag_coefs(jx, jy, jz, a, b, c, nullptr, nullptr, outloc);

          delft(jx, jz) = a * ft(jx - 1, jz) + b * ft(jx, jz) + c * ft(jx + 1, jz);
        }
      }

      // Reverse FFT
      for (int jx = localmesh->xstart; jx <= localmesh->xend; jx++) {

        irfft(&delft(jx, 0), ncz, &result(jx, jy, 0));
      }
    }
  } else {
    result = G1 * ::DDX(f, outloc) + G3 * ::DDZ(f, outloc) + g11 * ::D2DX2(f, outloc)
             + g33 * ::D2DZ2(f, outloc) + 2 * g13 * ::D2DXDZ(f, outloc);
  };

  ASSERT2(result.getLocation() == outloc);

  return result;
}

FieldPerp Coordinates::Delp2(const FieldPerp& f, CELL_LOC outloc, bool useFFT) {
  TRACE("Coordinates::Delp2( FieldPerp )");

  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }

  ASSERT1(location == outloc);
  ASSERT1(f.getLocation() == outloc);

  if (localmesh->GlobalNx == 1 && localmesh->GlobalNz == 1) {
    // copy mesh, location, etc
    return f * 0;
  }
  ASSERT2(localmesh->xstart > 0); // Need at least one guard cell

  FieldPerp result{emptyFrom(f).setLocation(outloc)};

  int jy = f.getIndex();
  result.setIndex(jy);

  if (useFFT) {
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
      for (int jx = localmesh->xstart; jx <= localmesh->xend; jx++) {
        // Perform x derivative

        dcomplex a, b, c;
        laplace_tridag_coefs(jx, jy, jz, a, b, c);

        delft(jx, jz) = a * ft(jx - 1, jz) + b * ft(jx, jz) + c * ft(jx + 1, jz);
      }
    }

    // Reverse FFT
    for (int jx = localmesh->xstart; jx <= localmesh->xend; jx++) {
      irfft(&delft(jx, 0), ncz, &result(jx, 0));
    }

  } else {
    throw BoutException("Non-fourier Delp2 not currently implented for FieldPerp.");
    // Would be the following but don't have standard derivative operators for FieldPerps
    // yet
    // result = G1 * ::DDX(f, outloc) + G3 * ::DDZ(f, outloc) + g11 * ::D2DX2(f, outloc)
    //          + g33 * ::D2DZ2(f, outloc) + 2 * g13 * ::D2DXDZ(f, outloc);
  };

  return result;
}

Field2D Coordinates::Laplace_par(const Field2D& f, CELL_LOC outloc) {
  ASSERT1(location == outloc || outloc == CELL_DEFAULT);
  return D2DY2(f, outloc) / g_22 + DDY(J / g_22, outloc) * DDY(f, outloc) / J;
}

Field3D Coordinates::Laplace_par(const Field3D& f, CELL_LOC outloc) {
  ASSERT1(location == outloc || outloc == CELL_DEFAULT);
  return D2DY2(f, outloc) / g_22 + DDY(J / g_22, outloc) * ::DDY(f, outloc) / J;
}

// Full Laplacian operator on scalar field

Field2D Coordinates::Laplace(const Field2D& f, CELL_LOC outloc) {
  TRACE("Coordinates::Laplace( Field2D )");
  ASSERT1(location == outloc || outloc == CELL_DEFAULT);

  Field2D result = G1 * DDX(f, outloc) + G2 * DDY(f, outloc) + g11 * D2DX2(f, outloc)
                   + g22 * D2DY2(f, outloc) + 2.0 * g12 * D2DXDY(f, outloc);

  return result;
}

Field3D Coordinates::Laplace(const Field3D& f, CELL_LOC outloc) {
  TRACE("Coordinates::Laplace( Field3D )");
  ASSERT1(location == outloc || outloc == CELL_DEFAULT);

  Field3D result = G1 * ::DDX(f, outloc) + G2 * ::DDY(f, outloc) + G3 * ::DDZ(f, outloc)
                   + g11 * D2DX2(f, outloc) + g22 * D2DY2(f, outloc)
                   + g33 * D2DZ2(f, outloc)
                   + 2.0 * (g12 * D2DXDY(f, outloc) + g13 * D2DXDZ(f, outloc)
                            + g23 * D2DYDZ(f, outloc));

  return result;
}

// Full perpendicular Laplacian, in form of inverse of Laplacian operator in LaplaceXY
// solver
Field2D Coordinates::Laplace_perpXY(const Field2D& A, const Field2D& f) {
  TRACE("Coordinates::Laplace_perpXY( Field2D )");

  Field2D result;
  result.allocate();
  BoutReal thisA, thisJ, thisg11, thisdx, thisg_22, thisg23, thisg_23, thisdy, val;
  for (auto i : result.region(RGN_NOBNDRY)) {
    result[i] = 0.;

    // outer x boundary
    thisA = 0.5 * (A[i] + A[i.xp()]);
    thisJ = 0.5 * (J[i] + J[i.xp()]);
    thisg11 = 0.5 * (g11[i] + g11[i.xp()]);
    thisdx = 0.5 * (dx[i] + dx[i.xp()]);
    val = thisA * thisJ * thisg11 / (J[i] * thisdx * dx[i]);
    result[i] += val * (f[i.xp()] - f[i]);

    // inner x boundary
    thisA = 0.5 * (A[i] + A[i.xm()]);
    thisJ = 0.5 * (J[i] + J[i.xm()]);
    thisg11 = 0.5 * (g11[i] + g11[i.xm()]);
    thisdx = 0.5 * (dx[i] + dx[i.xm()]);
    val = thisA * thisJ * thisg11 / (J[i] * thisdx * dx[i]);
    result[i] += val * (f[i.xm()] - f[i]);

    // upper y boundary
    thisA = 0.5 * (A[i] + A[i.yp()]);
    thisJ = 0.5 * (J[i] + J[i.yp()]);
    thisg_22 = 0.5 * (g_22[i] + g_22[i.yp()]);
    thisg23 = 0.5 * (g23[i] + g23[i.yp()]);
    thisg_23 = 0.5 * (g_23[i] + g_23[i.yp()]);
    thisdy = 0.5 * (dy[i] + dy[i.yp()]);
    val = -thisA * thisJ * thisg23 * thisg_23 / (thisg_22 * J[i] * thisdy * dy[i]);
    result[i] += val * (f[i.yp()] - f[i]);

    // lower y boundary
    thisA = 0.5 * (A[i] + A[i.ym()]);
    thisJ = 0.5 * (J[i] + J[i.ym()]);
    thisg_22 = 0.5 * (g_22[i] + g_22[i.ym()]);
    thisg23 = 0.5 * (g23[i] + g23[i.ym()]);
    thisg_23 = 0.5 * (g_23[i] + g_23[i.ym()]);
    thisdy = 0.5 * (dy[i] + dy[i.ym()]);
    val = -thisA * thisJ * thisg23 * thisg_23 / (thisg_22 * J[i] * thisdy * dy[i]);
    result[i] += val * (f[i.ym()] - f[i]);
  }

  return result;
}
