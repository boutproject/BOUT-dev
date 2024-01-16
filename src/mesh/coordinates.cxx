/**************************************************************************
 * Differential geometry
 * Calculates the covariant metric tensor, and christoffel symbol terms
 * given the contravariant metric tensor terms
 **************************************************************************/

#include <bout/assert.hxx>
#include <bout/constants.hxx>
#include <bout/coordinates.hxx>
#include <bout/msg_stack.hxx>
#include <bout/sys/timer.hxx>
#include <bout/utils.hxx>

#include <bout/fft.hxx>

#include "parallel/fci.hxx"
#include "parallel/shiftedmetricinterp.hxx"

#include <bout/bout_types.hxx>

// use anonymous namespace so this utility function is not available outside this file
namespace {

template <typename T, typename... Ts>
// Use sendY()/sendX() and wait() instead of Mesh::communicate() to ensure we
// don't try to calculate parallel slices as Coordinates are not constructed yet
void communicate(T& t, Ts... ts) {
  FieldGroup g(t, ts...);
  auto h = t.getMesh()->sendY(g);
  t.getMesh()->wait(h);
  h = t.getMesh()->sendX(g);
  t.getMesh()->wait(h);
}

#if BOUT_USE_METRIC_3D
Field3D interpolateAndExtrapolate(const Field3D& f_, CELL_LOC location,
                                  bool extrapolate_x, bool extrapolate_y,
                                  bool no_extra_interpolate, ParallelTransform* pt_) {

  Mesh* localmesh = f_.getMesh();
  Field3D result;
  Field3D f = f_;
  ParallelTransform* pt_f;
  if (f.getCoordinates() == nullptr) {
    // if input f is member of the Coordinates we are currently constructing, it will not
    // have Coordinates and needs to use the passed-in ParallelTransform
    pt_f = pt_;
  } else {
    // if input f is from Coordinates at a different location, it will have its own
    // Coordinates, and we should use its ParallelTransform
    pt_f = &f.getCoordinates()->getParallelTransform();
  }
  if (f.getDirectionY() != YDirectionType::Standard) {
    if (pt_f->canToFromFieldAligned()) {
      f = pt_f->fromFieldAligned(f);
    } else {
      f.setDirectionY(YDirectionType::Standard);
    }
  }
  if (location == CELL_YLOW and f.getLocation() != CELL_YLOW) {
    auto f_aligned = pt_f->toFieldAligned(f, "RGN_NOX");
    result = interp_to(f_aligned, location, "RGN_NOBNDRY");
    ParallelTransform* pt_result;
    if (result.getCoordinates() == nullptr) {
      pt_result = pt_;
    } else {
      pt_result = &result.getCoordinates()->getParallelTransform();
    }
    result = pt_result->fromFieldAligned(result, "RGN_NOBNDRY");
  } else {
    result = interp_to(f, location, "RGN_NOBNDRY");
  }
  // Ensure result's data is unique. Otherwise result might be a duplicate of
  // f (if no interpolation is needed, e.g. if interpolation is in the
  // z-direction); then f would be communicated. Since this function is used
  // on geometrical quantities that might not be periodic in y even on closed
  // field lines (due to dependence on integrated shear), we don't want to
  // communicate f. We will sort out result's boundary guard cells below, but
  // not f's so we don't want to change f.
  result.allocate();
  communicate(result);

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
          for (int zi = 0; zi < localmesh->LocalNz; ++zi) {
            result(bndry->x, bndry->y, zi) =
                (9.
                     * (f(bndry->x - bndry->bx, bndry->y - bndry->by, zi)
                        + f(bndry->x, bndry->y, zi))
                 - f(bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, zi)
                 - f(bndry->x + bndry->bx, bndry->y + bndry->by, zi))
                / 16.;
          }
        }
        // set boundary guard cells
        if ((bndry->bx != 0 && localmesh->GlobalNx - 2 * bndry->width >= 3)
            || (bndry->by != 0 && localmesh->GlobalNy - 2 * bndry->width >= 3)) {
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
            for (int zi = 0; zi < localmesh->LocalNz; ++zi) {
              result(xi, yi, zi) =
                  3.0 * result(xi - bndry->bx, yi - bndry->by, zi)
                  - 3.0 * result(xi - 2 * bndry->bx, yi - 2 * bndry->by, zi)
                  + result(xi - 3 * bndry->bx, yi - 3 * bndry->by, zi);
            }
          }
        } else {
          // not enough grid points to extrapolate, set equal to last grid point
          for (int i = extrap_start; i < bndry->width; i++) {
            for (int zi = 0; zi < localmesh->LocalNz; ++zi) {
              result(bndry->x + i * bndry->bx, bndry->y + i * bndry->by, zi) =
                  result(bndry->x - bndry->bx, bndry->y - bndry->by, zi);
            }
          }
        }
      }
    }
  }
#if CHECK > 0
  if (not(
          // if include_corner_cells=true, then we extrapolate valid data into the
          // corner cells if they are not already filled
          localmesh->include_corner_cells

          // if we are not extrapolating at all, the corner cells should contain valid
          // data
          or (not extrapolate_x and not extrapolate_y))) {
    // Invalidate corner guard cells
    for (int i = 0; i < localmesh->xstart; i++) {
      for (int j = 0; j < localmesh->ystart; j++) {
        for (int k = 0; k < localmesh->LocalNz; ++k) {
          result(i, j, k) = BoutNaN;
          result(i, localmesh->LocalNy - 1 - j, k) = BoutNaN;
          result(localmesh->LocalNx - 1 - i, j, k) = BoutNaN;
          result(localmesh->LocalNx - 1 - i, localmesh->LocalNy - 1 - j, k) = BoutNaN;
        }
      }
    }
  }
#endif // CHECK > 0

  return result;
}
#endif // BOUT_USE_METRIC_3D

// If the CELL_CENTRE variable was read, the staggered version is required to
// also exist for consistency
void checkStaggeredGet(Mesh* mesh, const std::string& name, const std::string& suffix) {
  if (mesh->sourceHasVar(name) != mesh->sourceHasVar(name + suffix)) {
    throw BoutException("Attempting to read staggered fields from grid, but " + name
                        + " is not present in both CELL_CENTRE and staggered versions.");
  }
}

// convenience function for repeated code
int getAtLoc(Mesh* mesh, Coordinates::FieldMetric& var, const std::string& name,
             const std::string& suffix, CELL_LOC location, BoutReal default_value = 0.) {

  checkStaggeredGet(mesh, name, suffix);
  return mesh->get(var, name + suffix, default_value, false, location);
}

auto getAtLoc(Mesh* mesh, const std::string& name, const std::string& suffix,
              CELL_LOC location, BoutReal default_value = 0.) {

  checkStaggeredGet(mesh, name, suffix);
  return mesh->get(name + suffix, default_value, false, location);
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
    // in 2D metric, same as CELL_CENTRE
    return bout::build::use_metric_3d ? "_zlow" : "";
  }
  default: {
    throw BoutException(
        "Incorrect location passed to "
        "Coordinates(Mesh*,const CELL_LOC,const Coordinates*) constructor.");
  }
  }
}

} // anonymous namespace

Coordinates::FieldMetric Coordinates::getAtLocOrUnaligned(Mesh* mesh,
                                                          const std::string& name,
                                                          BoutReal default_value,
                                                          const std::string& suffix,
                                                          CELL_LOC cell_location) {
  if (cell_location == CELL_CENTRE) {
    return getUnaligned(name, default_value);
  }
  // grid data source has staggered fields, so read instead of interpolating
  // Diagonal components of metric tensor g^{ij} (default to 1)
  return getAtLoc(mesh, name, suffix, cell_location, default_value);
}

Coordinates::FieldMetric Coordinates::getUnaligned(const std::string& name,
                                                   BoutReal default_value) {

  auto field = localmesh->get(name, default_value, false);
  if (field.getDirectionY() == YDirectionType::Aligned
      and transform->canToFromFieldAligned()) {
    return transform->fromFieldAligned(field);
  }
  field.setDirectionY(YDirectionType::Standard);
  return field;
}

Coordinates::FieldMetric Coordinates::getUnalignedAtLocationAndFillGuards(
    Mesh* mesh, const std::string& name, BoutReal default_value,
    const std::string& suffix, CELL_LOC cell_location, bool extrapolate_x,
    bool extrapolate_y, bool no_extra_interpolate,
    ParallelTransform* pParallelTransform) {

  auto field = getAtLocOrUnaligned(mesh, name, default_value, suffix, cell_location);
  if (suffix.empty()) {
    no_extra_interpolate = false;
    pParallelTransform = transform.get();
  }

  return field.getMesh()->interpolateAndExtrapolate(field, cell_location, extrapolate_x,
                                                    extrapolate_y, no_extra_interpolate,
                                                    pParallelTransform);
}

Coordinates::Coordinates(Mesh* mesh, FieldMetric dx, FieldMetric dy, FieldMetric dz,
                         FieldMetric J, FieldMetric Bxy, const FieldMetric& g11,
                         const FieldMetric& g22, const FieldMetric& g33,
                         const FieldMetric& g12, const FieldMetric& g13,
                         const FieldMetric& g23, const FieldMetric& g_11,
                         const FieldMetric& g_22, const FieldMetric& g_33,
                         const FieldMetric& g_12, const FieldMetric& g_13,
                         const FieldMetric& g_23, FieldMetric ShiftTorsion,
                         FieldMetric IntShiftTorsion)
    : nz(mesh->LocalNz), localmesh(mesh), location(CELL_CENTRE),
      differential_operators(mesh->getDifferentialOperators()), dx_(std::move(dx)),
      dy_(std::move(dy)), dz_(std::move(dz)), ShiftTorsion_(std::move(ShiftTorsion)),
      IntShiftTorsion_(std::move(IntShiftTorsion)),
      contravariantMetricTensor(g11, g22, g33, g12, g13, g23),
      covariantMetricTensor(g_11, g_22, g_33, g_12, g_13, g_23), J_(std::move(J)),
      Bxy_(std::move(Bxy)){ASSERT0(differential_operators != nullptr)};

Coordinates::Coordinates(Mesh* mesh, Options* options, const CELL_LOC loc,
                         const Coordinates* coords_in, bool force_interpolate_from_centre)
    : localmesh(mesh), location(loc),
      differential_operators(mesh->getDifferentialOperators()), dx_(1., mesh),
      dy_(1., mesh), dz_(1., mesh), d1_dx_(mesh), d1_dy_(mesh), d1_dz_(mesh),
      ShiftTorsion_(mesh), IntShiftTorsion_(mesh),
      contravariantMetricTensor(1., 1., 1., 0, 0, 0, mesh),
      // Identity metric tensor
      covariantMetricTensor(1., 1., 1., 0, 0, 0, mesh), J_(1., mesh), Bxy_(1., mesh) {
  ASSERT0(differential_operators != nullptr);

  if (options == nullptr) {
    options = Options::getRoot()->getSection("mesh");
  }

  nz = mesh->LocalNz;

  const std::string suffix = getLocationSuffix(location);

  if (coords_in == nullptr || suffix.empty()
      || (!force_interpolate_from_centre && mesh->sourceHasVar("dx" + suffix))) {
    setBoundaryCells(mesh, options, coords_in, suffix);
  } else {
    interpolateFieldsFromOtherCoordinates(mesh, options, coords_in);
  }
}

void Coordinates::interpolateFieldsFromOtherCoordinates(const Mesh* mesh,
                                                        Options* options,
                                                        const Coordinates* coords_in) {

  std::function<const FieldMetric(const FieldMetric)> const
      interpolateAndExtrapolate_function = [this](const FieldMetric& component) {
        return localmesh->interpolateAndExtrapolate(component, location, true, true,
                                                    false, transform.get());
      };

  const auto region = std::basic_string("RGN_NOBNDRY");
  setContravariantMetricTensor(coords_in->getContravariantMetricTensor(), region);

  applyToContravariantMetricTensor(interpolateAndExtrapolate_function);
  applyToCovariantMetricTensor(interpolateAndExtrapolate_function);

  // Check input metrics
  checkContravariant();
  checkCovariant();

  setJ(localmesh->interpolateAndExtrapolate(coords_in->J(), location, true, true, false,
                                            transform.get()));
  setBxy(localmesh->interpolateAndExtrapolate(coords_in->Bxy(), location, true, true,
                                              false, transform.get()));

  bout::checkFinite(J(), "The Jacobian", "RGN_NOCORNERS");
  bout::checkPositive(J(), "The Jacobian", "RGN_NOCORNERS");
  bout::checkFinite(Bxy(), "Bxy", "RGN_NOCORNERS");
  bout::checkPositive(Bxy(), "Bxy", "RGN_NOCORNERS");

  ShiftTorsion_ = localmesh->interpolateAndExtrapolate(
      coords_in->ShiftTorsion(), location, true, true, false, transform.get());

  if (mesh->IncIntShear) {
    IntShiftTorsion_ = localmesh->interpolateAndExtrapolate(
        coords_in->IntShiftTorsion(), location, true, true, false, transform.get());
  }

  if (isUniform(coords_in->dz())) {
    dz_ = coords_in->dz();
    dz_.setLocation(location);
  } else {
    throw BoutException("We are asked to transform dz to get dz before we "
                        "have a transform, which "
                        "might require dz!\nPlease provide a dz for the "
                        "staggered quantity!");
  }
  setParallelTransform(options);
  dx_ = localmesh->interpolateAndExtrapolate(coords_in->dx(), location, true, true, false,
                                             transform.get());
  dy_ = localmesh->interpolateAndExtrapolate(coords_in->dy(), location, true, true, false,
                                             transform.get());
  // not really needed - we have used dz already ...
  dz_ = localmesh->interpolateAndExtrapolate(coords_in->dz(), location, true, true, false,
                                             transform.get());
}

// Note: If boundary cells were not loaded from the grid file, use
// 'interpolateAndExtrapolate' to set them. Ensures that derivatives are
// smooth at all the boundaries.
void Coordinates::setBoundaryCells(Mesh* mesh, Options* options,
                                   const Coordinates* coords_in,
                                   const std::string& suffix) {

  const bool extrapolate_x =
      (*options)["extrapolate_x"].withDefault(not mesh->sourceHasXBoundaryGuards());
  const bool extrapolate_y =
      (*options)["extrapolate_y"].withDefault(not mesh->sourceHasYBoundaryGuards());

  if (extrapolate_x) {
    output_warn.write(_("WARNING: extrapolating input mesh quantities into x-boundary "
                        "cells. Set option extrapolate_x=false to disable this.\n"));
  }

  if (extrapolate_y) {
    output_warn.write(_("WARNING: extrapolating input mesh quantities into y-boundary "
                        "cells. Set option extrapolate_y=false to disable this.\n"));
  }

  if (coords_in == nullptr) {
    mesh->get(dx_, "dx", 1.0, false);
    mesh->get(dy_, "dy", 1.0, false);
  }

  nz = mesh->LocalNz;

  {
    auto& options = Options::root();
    const bool has_zperiod = options.isSet("zperiod");
    const auto zmin = has_zperiod ? 0.0 : options["ZMIN"].withDefault(0.0);
    const auto zmax = has_zperiod ? 1.0 / options["zperiod"].withDefault(1.0)
                                  : options["ZMAX"].withDefault(1.0);

    const auto default_dz = (zmax - zmin) * TWOPI / nz;
    getAtLoc(mesh, dz_, "dz", suffix, location, default_dz);
  }

  // required early for differentiation.
  setParallelTransform(options);

  dz_ = localmesh->interpolateAndExtrapolate(dz_, location, extrapolate_x, extrapolate_y,
                                             false, transform.get());

  dx_ = getUnalignedAtLocationAndFillGuards(mesh, "dx", 1.0, suffix, location,
                                            extrapolate_x, extrapolate_y, false,
                                            transform.get());

  if (mesh->periodicX) {
    communicate(dx_);
  }

  dy_ = getUnalignedAtLocationAndFillGuards(mesh, "dy", 1.0, suffix, location,
                                            extrapolate_x, extrapolate_y, false,
                                            transform.get());

  // grid data source has staggered fields, so read instead of interpolating
  // Diagonal components of metric tensor g^{ij} (default to 1)
  //    TODO: Method `getAtLocAndFillGuards` violates command–query separation principle?
  FieldMetric g11, g22, g33, g12, g13, g23;

  // Diagonal components of metric tensor g^{ij} (default to 1)
  g11 = getUnalignedAtLocationAndFillGuards(mesh, "g11", 1.0, suffix, location,
                                            extrapolate_x, extrapolate_y, false,
                                            transform.get());
  g22 = getUnalignedAtLocationAndFillGuards(mesh, "g22", 1.0, suffix, location,
                                            extrapolate_x, extrapolate_y, false,
                                            transform.get());
  g33 = getUnalignedAtLocationAndFillGuards(mesh, "g33", 1.0, suffix, location,
                                            extrapolate_x, extrapolate_y, false,
                                            transform.get());
  // Off-diagonal elements. Default to 0
  g12 = getUnalignedAtLocationAndFillGuards(mesh, "g12", 0.0, suffix, location,
                                            extrapolate_x, extrapolate_y, false,
                                            transform.get());
  g13 = getUnalignedAtLocationAndFillGuards(mesh, "g13", 0.0, suffix, location,
                                            extrapolate_x, extrapolate_y, false,
                                            transform.get());
  g23 = getUnalignedAtLocationAndFillGuards(mesh, "g23", 0.0, suffix, location,
                                            extrapolate_x, extrapolate_y, false,
                                            transform.get());

  contravariantMetricTensor.setMetricTensor(MetricTensor(g11, g22, g33, g12, g13, g23));

  // Check input metrics
  checkContravariant();

  /// Find covariant metric components
  auto covariant_component_names = {"g_11", "g_22", "g_33", "g_12", "g_13", "g_23"};
  auto source_has_component = [&suffix, &mesh](const std::string& name) {
    return mesh->sourceHasVar(name + suffix);
  };

  // Check that all components are present
  if (std::all_of(begin(covariant_component_names), end(covariant_component_names),
                  source_has_component)) {

    FieldMetric g_11, g_22, g_33, g_12, g_13, g_23;
    g_11 = getAtLocOrUnaligned(mesh, "g_11", 1.0, suffix, location);
    g_22 = getAtLocOrUnaligned(mesh, "g_22", 1.0, suffix, location);
    g_33 = getAtLocOrUnaligned(mesh, "g_33", 1.0, suffix, location);
    g_12 = getAtLocOrUnaligned(mesh, "g_12", 0.0, suffix, location);
    g_13 = getAtLocOrUnaligned(mesh, "g_13", 0.0, suffix, location);
    g_23 = getAtLocOrUnaligned(mesh, "g_23", 0.0, suffix, location);
    covariantMetricTensor.setMetricTensor(
        MetricTensor(g_11, g_22, g_33, g_12, g_13, g_23));

    output_warn.write("\tWARNING! Covariant components of metric tensor set manually. "
                      "Contravariant components NOT recalculated\n");
  } else {
    covariantMetricTensor.setMetricTensor(contravariantMetricTensor.inverse());
    output_warn.write("Not all covariant components of metric tensor found. "
                      "Calculating all from the contravariant tensor\n");
  }

  // More robust to extrapolate derived quantities directly, rather than
  // deriving from extrapolated covariant metric components

  std::function<const FieldMetric(const FieldMetric)> const
      interpolateAndExtrapolate_function = [this, extrapolate_y,
                                            extrapolate_x](const FieldMetric& component) {
        return localmesh->interpolateAndExtrapolate(
            component, location, extrapolate_x, extrapolate_y, false, transform.get());
      };
  applyToCovariantMetricTensor(interpolateAndExtrapolate_function);

  // Check covariant metrics
  checkCovariant();

  /// Calculate Jacobian and Bxy
  jacobian();

  // Attempt to read J from the grid file
  if (!localmesh->sourceHasVar("J" + suffix)) {
    output_warn.write(
        "\tWARNING: Jacobian 'J_{:s}' not found. Calculating from metric tensor\n",
        suffix);
  } else {
    const auto Jcalc = getAtLoc(mesh, "J", suffix, location);
    setJ(localmesh->interpolateAndExtrapolate(Jcalc, location, extrapolate_x,
                                              extrapolate_y, false, transform.get()));

    // Compare calculated and loaded values
    output_warn.write("\tMaximum difference in J is {:e}\n", max(abs(J() - Jcalc)));

    auto J_value = J(); // TODO: There may be a better way
    communicate(J_value);

    // Re-evaluate Bxy using new J
    setBxy(sqrt(g_22()) / J());
  }

  // Check jacobian
  bout::checkFinite(J(), "J" + suffix, "RGN_NOCORNERS");
  bout::checkPositive(J(), "J" + suffix, "RGN_NOCORNERS");
  if (min(abs(J())) < 1.0e-10) {
    throw BoutException("\tERROR: Jacobian{:s} becomes very small\n", suffix);
  }

  // Attempt to read Bxy from the grid file
  if (!localmesh->sourceHasVar("Bxy" + suffix)) {
    output_warn.write("\tWARNING: Magnitude of B field 'Bxy_{:s}' not found. "
                      "Calculating from metric tensor\n",
                      suffix);
  } else {
    const auto Bcalc = getAtLoc(mesh, "Bxy", suffix, location);
    setBxy(localmesh->interpolateAndExtrapolate(Bcalc, location, extrapolate_x,
                                                extrapolate_y, false, transform.get()));
    output_warn.write("\tMaximum difference in Bxy is {:e}\n", max(abs(Bxy() - Bcalc)));
  }
  // Check Bxy
  bout::checkFinite(Bxy(), "Bxy" + suffix, "RGN_NOCORNERS");
  bout::checkPositive(Bxy(), "Bxy" + suffix, "RGN_NOCORNERS");

  if (getAtLoc(mesh, ShiftTorsion_, "ShiftTorsion", suffix, location) != 0) {
    output_warn.write("\tWARNING: No Torsion specified for zShift. "
                      "Derivatives may not be correct\n");
    ShiftTorsion_ = 0.0;
  }
  ShiftTorsion_ = localmesh->interpolateAndExtrapolate(
      ShiftTorsion_, location, extrapolate_x, extrapolate_y, false, transform.get());

  //////////////////////////////////////////////////////

  if (mesh->IncIntShear) {
    checkStaggeredGet(mesh, "IntShiftTorsion", suffix);
    if (mesh->get(IntShiftTorsion_, "IntShiftTorsion" + suffix, 0.0, false) != 0) {
      output_warn.write("\tWARNING: No Integrated torsion specified\n");
      IntShiftTorsion_ = 0.0;
    }
    IntShiftTorsion_.setLocation(location);
    IntShiftTorsion_ = localmesh->interpolateAndExtrapolate(
        IntShiftTorsion_, location, extrapolate_x, extrapolate_y, false, transform.get());
  } else {
    // IntShiftTorsion will not be used, but set to zero to avoid uninitialized field
    IntShiftTorsion_ = 0.;
  }
}

void Coordinates::outputVars(Options& output_options) {
  Timer const time("io");
  const std::string loc_string =
      (location == CELL_CENTRE) ? "" : "_" + toString(location);

  output_options["dx" + loc_string].force(dx(), "Coordinates");
  output_options["dy" + loc_string].force(dy(), "Coordinates");
  output_options["dz" + loc_string].force(dz(), "Coordinates");

  output_options["g11" + loc_string].force(g11(), "Coordinates");
  output_options["g22" + loc_string].force(g22(), "Coordinates");
  output_options["g33" + loc_string].force(g33(), "Coordinates");
  output_options["g12" + loc_string].force(g12(), "Coordinates");
  output_options["g13" + loc_string].force(g13(), "Coordinates");
  output_options["g23" + loc_string].force(g23(), "Coordinates");

  output_options["g_11" + loc_string].force(g_11(), "Coordinates");
  output_options["g_22" + loc_string].force(g_22(), "Coordinates");
  output_options["g_33" + loc_string].force(g_33(), "Coordinates");
  output_options["g_12" + loc_string].force(g_12(), "Coordinates");
  output_options["g_13" + loc_string].force(g_13(), "Coordinates");
  output_options["g_23" + loc_string].force(g_23(), "Coordinates");

  output_options["J" + loc_string].force(J(), "Coordinates");
  output_options["Bxy" + loc_string].force(Bxy(), "Coordinates");

  output_options["G1" + loc_string].force(G1(), "Coordinates");
  output_options["G2" + loc_string].force(G2(), "Coordinates");
  output_options["G3" + loc_string].force(G3(), "Coordinates");

  getParallelTransform().outputVars(output_options);
}

const Field2D& Coordinates::zlength() const {
  BOUT_OMP(critical)
  if (not zlength_cache) {
    zlength_cache = std::make_unique<Field2D>(0., localmesh);

#if BOUT_USE_METRIC_3D
    BOUT_FOR_SERIAL(i, dz.getRegion("RGN_ALL")) { (*zlength_cache)[i] += dz[i]; }
#else
    (*zlength_cache) = dz_ * nz;
#endif
  }

  return *zlength_cache;
}

int Coordinates::communicateAndCheckMeshSpacing() const {
  TRACE("Coordinates::communicateAndCheckMeshSpacing");

  auto tmp = dx(); // TODO: There must be a better way than this!
  communicate(tmp, dy(), dz(), g11(), g22(), g33(), g12(), g13(), g23(), g_11(), g_22(),
              g_33(), g_12(), g_13(), g_23(), J(), Bxy());

  output_progress.write("Calculating differential geometry terms\n");

  if (min(abs(dx())) < 1e-8) {
    throw BoutException("dx magnitude less than 1e-8");
  }

  if (min(abs(dy())) < 1e-8) {
    throw BoutException("dy magnitude less than 1e-8");
  }

  if (min(abs(dz())) < 1e-8) {
    throw BoutException("dz magnitude less than 1e-8");
  }

  return 0;
}

void Coordinates::recalculateAndReset(bool recalculate_staggered,
                                      bool force_interpolate_from_centre) {

  // Check input metrics
  checkContravariant();
  checkCovariant();

  calculateCommunicateAndExtrapolateChristoffelSymbols();

  correctionForNonUniformMeshes(force_interpolate_from_centre);

  if (location == CELL_CENTRE && recalculate_staggered) {
    // Re-calculate interpolated Coordinates at staggered locations
    localmesh->recalculateStaggeredCoordinates();
  }

  invalidateAndRecalculateCachedVariables();
}

void Coordinates::correctionForNonUniformMeshes(bool force_interpolate_from_centre) {
  OPTION(Options::getRoot(), non_uniform_, true);

  FieldMetric d2x(localmesh);
  FieldMetric d2y(localmesh);
  FieldMetric const d2z(localmesh); // d^2 x / d i^2

  // Read correction for non-uniform meshes
  std::string const suffix = getLocationSuffix(location);

  auto extrapolate_x = true;
  auto extrapolate_y = true;
  if (location == CELL_CENTRE
      or (!force_interpolate_from_centre and localmesh->sourceHasVar("dx" + suffix))) {
    extrapolate_x = not localmesh->sourceHasXBoundaryGuards();
    extrapolate_y = not localmesh->sourceHasYBoundaryGuards();
  }

  if (localmesh->get(d2x, "d2x" + suffix, 0.0, false, location) != 0) {
    output_warn.write("\tWARNING: differencing quantity 'd2x' not found. "
                      "Calculating from dx\n");
    d1_dx_ = bout::derivatives::index::DDX(1. / dx()); // d/di(1/dx)

    communicate(d1_dx_);
    d1_dx_ = localmesh->interpolateAndExtrapolate(d1_dx_, location, true, true, true,
                                                  transform.get());
  } else {
    d2x.setLocation(location);
    // set boundary cells if necessary
    d2x = localmesh->interpolateAndExtrapolate(d2x, location, extrapolate_x,
                                               extrapolate_y, false, transform.get());

    d1_dx_ = -d2x / (dx() * dx());
  }

  if (localmesh->get(d2y, "d2y" + suffix, 0.0, false, location) != 0) {
    output_warn.write("\tWARNING: differencing quantity 'd2y' not found. "
                      "Calculating from dy\n");
    d1_dy_ = bout::derivatives::index::DDY(1. / dy()); // d/di(1/dy)

    communicate(d1_dy_);
    d1_dy_ = localmesh->interpolateAndExtrapolate(d1_dy_, location, true, true, true,
                                                  transform.get());
  } else {
    d2y.setLocation(location);
    // set boundary cells if necessary
    d2y = localmesh->interpolateAndExtrapolate(d2y, location, extrapolate_x,
                                               extrapolate_y, false, transform.get());

    d1_dy_ = -d2y / (dy() * dy());
  }

#if BOUT_USE_METRIC_3D
  if (localmesh->get(d2z, "d2z" + suffix, 0.0, false)) {
    output_warn.write("\tWARNING: differencing quantity 'd2z' not found. "
                      "Calculating from dz\n");
    d1_dz_ = bout::derivatives::index::DDZ(1. / dz());
    communicate(d1_dz_);
    d1_dz_ = localmesh->interpolateAndExtrapolate(d1_dz_, location, true, true, true,
                                                  transform.get());
  } else {
    d2z.setLocation(location);
    // set boundary cells if necessary
    d2z = localmesh->interpolateAndExtrapolate(d2z, location, extrapolate_x,
                                               extrapolate_y, false, transform.get());

    d1_dz_ = -d2z / (dz() * dz());
  }
#else
  d1_dz_ = 0;
#endif

  auto tmp = d1_dx(); // TODO: There must be a better way than this!
  communicate(tmp, d1_dy(), d1_dz());
}

void Coordinates::calculateCommunicateAndExtrapolateChristoffelSymbols() {

  christoffel_symbols_cache.reset();

  auto tmp = J() * g12();
  communicate(tmp);
  setG1((DDX(J() * g11()) + DDY(tmp) + DDZ(J() * g13())) / J());
  tmp = J() * g22();
  communicate(tmp);
  setG2((DDX(J() * g12()) + DDY(tmp) + DDZ(J() * g23())) / J());
  tmp = J() * g23();
  communicate(tmp);
  setG3((DDX(J() * g13()) + DDY(tmp) + DDZ(J() * g33())) / J());

  communicateChristoffelSymbolTerms();

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

  std::function<const FieldMetric(const FieldMetric)> const
      interpolateAndExtrapolate_function = [this](const FieldMetric& component) {
        return localmesh->interpolateAndExtrapolate(component, location, true, true,
                                                    false, transform.get());
      };

  applyToChristoffelSymbols(interpolateAndExtrapolate_function);
}

MetricTensor::FieldMetric Coordinates::recalculateJacobian() {

  // calculate Jacobian using g^-1 = det[g^ij], J = sqrt(g)
  auto g = contravariantMetricTensor.g11() * contravariantMetricTensor.g22()
               * contravariantMetricTensor.g33()
           + 2.0 * contravariantMetricTensor.g12() * contravariantMetricTensor.g13()
                 * contravariantMetricTensor.g23()
           - contravariantMetricTensor.g11() * contravariantMetricTensor.g23()
                 * contravariantMetricTensor.g23()
           - contravariantMetricTensor.g22() * contravariantMetricTensor.g13()
                 * contravariantMetricTensor.g13()
           - contravariantMetricTensor.g33() * contravariantMetricTensor.g12()
                 * contravariantMetricTensor.g12();

  // Check that g is positive
  bout::checkPositive(g, "The determinant of g^ij", "RGN_NOBNDRY");

  return 1. / sqrt(g);
}

MetricTensor::FieldMetric Coordinates::recalculateBxy() {

  return sqrt(covariantMetricTensor.g22()) / J_;
}

void Coordinates::jacobian() {
  TRACE("Coordinates::jacobian");

  const bool extrapolate_x = not localmesh->sourceHasXBoundaryGuards();
  const bool extrapolate_y = not localmesh->sourceHasYBoundaryGuards();

  try {

    const auto jacobian = recalculateJacobian();
    // More robust to extrapolate derived quantities directly, rather than
    // deriving from extrapolated covariant metric components
    setJ(localmesh->interpolateAndExtrapolate(jacobian, location, extrapolate_x,
                                              extrapolate_y, false));

    const auto Bxy = recalculateBxy();
    //    CELL_LOC location, ParallelTransform* pParallelTransform
    setBxy(localmesh->interpolateAndExtrapolate(Bxy, location, extrapolate_x,
                                                extrapolate_y, false, transform.get()));
  } catch (BoutException&) {
    output_error.write("\tError in jacobian call\n");
    throw;
  }
}

namespace {
// Utility function for fixing up guard cells of zShift
void fixZShiftGuards(Field2D& zShift) {
  auto* localmesh = zShift.getMesh();

  // extrapolate into boundary guard cells if necessary
  zShift = localmesh->interpolateAndExtrapolate(
      zShift, zShift.getLocation(), not localmesh->sourceHasXBoundaryGuards(),
      not localmesh->sourceHasYBoundaryGuards(), false);

  // make sure zShift has been communicated
  communicate(zShift);

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
} // namespace

void Coordinates::setParallelTransform(Options* options) {
  auto* ptoptions = options->getSection("paralleltransform");

  std::string ptstr;
  ptoptions->get("type", ptstr, "identity");

  // Convert to lower case for comparison
  ptstr = lowercase(ptstr);

  if (ptstr == "identity") {
    // Identity method i.e. no transform needed
    transform =
        bout::utils::make_unique<ParallelTransformIdentity>(*localmesh, ptoptions);

  } else if (ptstr == "shifted" or ptstr == "shiftedinterp") {
    // Shifted metric method

    Field2D zShift{localmesh};

    // Read the zShift angle from the mesh
    std::string const suffix = getLocationSuffix(location);
    if (localmesh->sourceHasVar("dx" + suffix)) {
      // Grid file has variables at this location, so should be able to read
      checkStaggeredGet(localmesh, "zShift", suffix);
      if (localmesh->get(zShift, "zShift" + suffix, 0.0, false, location) != 0) {
        // No zShift variable. Try qinty in BOUT grid files
        if (localmesh->get(zShift, "qinty" + suffix, 0.0, false, location) != 0) {
          // Failed to find either variable, cannot use ShiftedMetric
          throw BoutException("Could not read zShift" + suffix + " from grid file");
        }
      }
    } else {
      if (location == CELL_YLOW and bout::build::use_metric_3d) {
        throw BoutException("Cannot interpolate zShift to construct ShiftedMetric when "
                            "using 3d metrics. You must provide zShift_ylow in the grid "
                            "file.");
      }
      Field2D zShift_centre;
      if (localmesh->get(zShift_centre, "zShift", 0.0, false) != 0) {
        // No zShift variable. Try qinty in BOUT grid files
        if (localmesh->get(zShift_centre, "qinty", 0.0, false) != 0) {
          // Failed to find either variable, cannot use ShiftedMetric
          throw BoutException("Could not read zShift from grid file");
        }
      }

      fixZShiftGuards(zShift_centre);

      zShift = localmesh->interpolateAndExtrapolate(zShift_centre, location, true, true,
                                                    false, transform.get());
    }

    fixZShiftGuards(zShift);

    if (ptstr == "shifted") {
      transform = bout::utils::make_unique<ShiftedMetric>(*localmesh, location, zShift,
                                                          getUniform(zlength()));
    } else if (ptstr == "shiftedinterp") {
      transform = bout::utils::make_unique<ShiftedMetricInterp>(
          *localmesh, location, zShift, getUniform(zlength()));
    }

  } else if (ptstr == "fci") {

    if (location != CELL_CENTRE) {
      throw BoutException("FCITransform is not available on staggered grids.");
    }

    // Flux Coordinate Independent method
    const bool fci_zperiodic = (*ptoptions)["z_periodic"].withDefault(true);
    transform = bout::utils::make_unique<FCITransform>(*localmesh, dy(), fci_zperiodic,
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
                         const std::string& region) const {
  ASSERT1(location == loc || loc == CELL_DEFAULT)
  return differential_operators->DDX(f, dx(), loc, method, region);
}

Field3D Coordinates::DDX(const Field3D& f, CELL_LOC outloc, const std::string& method,
                         const std::string& region) const {
  return differential_operators->DDX(f, dx(), dz(), IntShiftTorsion_, outloc, method,
                                     region);
}

Field2D Coordinates::DDY(const Field2D& f, CELL_LOC loc, const std::string& method,
                         const std::string& region) const {
  ASSERT1(location == loc || loc == CELL_DEFAULT)
  return differential_operators->DDY(f, dy(), loc, method, region);
}

Field3D Coordinates::DDY(const Field3D& f, CELL_LOC outloc, const std::string& method,
                         const std::string& region) const {
  return differential_operators->DDY(f, dy(), outloc, method, region);
}

Field2D Coordinates::DDZ(const Field2D& f, CELL_LOC loc, const std::string& method,
                         const std::string& region) const {
  ASSERT1(location == loc || loc == CELL_DEFAULT)
  ASSERT1(f.getMesh() == localmesh)
  if (loc == CELL_DEFAULT) {
    loc = f.getLocation();
  }
  return differential_operators->DDZ(f, loc, method, region);
}

Field3D Coordinates::DDZ(const Field3D& f, CELL_LOC outloc, const std::string& method,
                         const std::string& region) const {
  return differential_operators->DDZ(f, dz(), outloc, method, region);
}

/////////////////////////////////////////////////////////
// Parallel gradient

Field2D Coordinates::Grad_par(const Field2D& var, MAYBE_UNUSED(CELL_LOC outloc),
                              const std::string& UNUSED(method)) {
  TRACE("Coordinates::Grad_par( Field2D )");

  return differential_operators->Grad_par(var, dy());
}

Field3D Coordinates::Grad_par(const Field3D& var, CELL_LOC outloc,
                              const std::string& method) {
  TRACE("Coordinates::Grad_par( Field3D )");
  ASSERT1(location == outloc || outloc == CELL_DEFAULT)

  return differential_operators->Grad_par(var, dy(), outloc, method);
}

/////////////////////////////////////////////////////////
// Vpar_Grad_par
// vparallel times the parallel derivative along unperturbed B-field

Field2D Coordinates::Vpar_Grad_par(const Field2D& v, const Field2D& f,
                                   MAYBE_UNUSED(CELL_LOC outloc),
                                   const std::string& UNUSED(method)) {
  ASSERT1(location == outloc || (outloc == CELL_DEFAULT && location == f.getLocation()))

  return differential_operators->Vpar_Grad_par(v, f, this, outloc);
}

Field3D Coordinates::Vpar_Grad_par(const Field3D& v, const Field3D& f, CELL_LOC outloc,
                                   const std::string& method) {
  ASSERT1(location == outloc || outloc == CELL_DEFAULT)

  return differential_operators->Vpar_Grad_par(v, f, this, outloc, method);
}

/////////////////////////////////////////////////////////
// Parallel divergence

Field2D Coordinates::Div_par(const Field2D& f, CELL_LOC outloc,
                             const std::string& method) {
  TRACE("Coordinates::Div_par( Field2D )");
  ASSERT1(location == outloc || outloc == CELL_DEFAULT)

  return differential_operators->Div_par(f, Bxy(), dy(), outloc, method);
}

Field3D Coordinates::Div_par(const Field3D& f, CELL_LOC outloc,
                             const std::string& method) {
  TRACE("Coordinates::Div_par( Field3D )");
  ASSERT1(location == outloc || outloc == CELL_DEFAULT)

  // Need Bxy at location of f, which might be different from location of this
  // Coordinates object
  auto Bxy_floc = f.getCoordinates()->Bxy();

  if (!f.hasParallelSlices()) {
    // No yup/ydown fields. The Grad_par operator will
    // shift to field aligned coordinates
    return differential_operators->Div_par(f, dy(), Bxy(), outloc, method);
  }

  // Need to modify yup and ydown fields
  Field3D f_B = f / Bxy_floc;
  f_B.splitParallelSlices();
  for (int i = 0; i < f.getMesh()->ystart; ++i) {
    f_B.yup(i) = f.yup(i) / Bxy_floc.yup(i);
    f_B.ydown(i) = f.ydown(i) / Bxy_floc.ydown(i);
  }
  return Bxy() * Grad_par(f_B, outloc, method);
}

/////////////////////////////////////////////////////////
// second parallel derivative (b dot Grad)(b dot Grad)
// Note: For parallel Laplacian use Laplace_par

Field2D Coordinates::Grad2_par2(const Field2D& f, CELL_LOC outloc,
                                const std::string& method) {
  TRACE("Coordinates::Grad2_par2( Field2D )");
  ASSERT1(location == outloc || (outloc == CELL_DEFAULT && location == f.getLocation()))

  auto result = Grad2_par2_DDY_invSg(outloc, method) * DDY(f, outloc, method)
                + differential_operators->D2DY2(f, outloc, method) / g_22(); // derivs 209

  return result;
}

Field3D Coordinates::Grad2_par2(const Field3D& f, CELL_LOC outloc,
                                const std::string& method) {
  TRACE("Coordinates::Grad2_par2( Field3D )");
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  ASSERT1(location == outloc)

  Field3D result = differential_operators->DDY(f, outloc, method);

  Field3D r2 = differential_operators->D2DY2(f, outloc, method) / g_22();

  result = Grad2_par2_DDY_invSg(outloc, method) * result + r2;

  ASSERT2(result.getLocation() == outloc)

  return result;
}

/////////////////////////////////////////////////////////
// perpendicular Laplacian operator

#include <bout/invert_laplace.hxx> // Delp2 uses same coefficients as inversion code
#include <utility>

Field2D Coordinates::Delp2(const Field2D& f, CELL_LOC outloc, bool UNUSED(useFFT)) {
  TRACE("Coordinates::Delp2( Field2D )");
  ASSERT1(location == outloc || outloc == CELL_DEFAULT)

  return G1() * DDX(f, outloc) + g11() * differential_operators->D2DX2(f, outloc);
}

Field3D Coordinates::Delp2(const Field3D& f, CELL_LOC outloc, bool useFFT) {
  TRACE("Coordinates::Delp2( Field3D )");

  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }

  ASSERT1(location == outloc)
  ASSERT1(f.getLocation() == outloc)

  if (localmesh->GlobalNx == 1 && localmesh->GlobalNz == 1) {
    // copy mesh, location, etc
    return f * 0;
  }
  ASSERT2(localmesh->xstart > 0) // Need at least one guard cell

  Field3D result{emptyFrom(f).setLocation(outloc)};

  if (useFFT and not bout::build::use_metric_3d) {
    int ncz = localmesh->LocalNz;

    // Allocate memory
    auto ft = Matrix<dcomplex>(localmesh->LocalNx, ncz / 2 + 1);
    auto delft = Matrix<dcomplex>(localmesh->LocalNx, ncz / 2 + 1);

    // Loop over y indices
    // Note: should not include y-guard or y-boundary points here as that would
    // use values from corner cells in dx, which may not be initialised.
    for (int jy = localmesh->ystart; jy <= localmesh->yend; jy++) {

      // Take forward FFT

      for (int jx = 0; jx < localmesh->LocalNx; jx++) {
        rfft(&f(jx, jy, 0), ncz, &ft(jx, 0));
      }

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
    result = G1() * differential_operators->DDX(f, outloc)
             + G3() * differential_operators->DDZ(f, outloc)
             + g11() * differential_operators->D2DX2(f, outloc)
             + g33() * differential_operators->D2DZ2(f, outloc)
             + 2 * g13() * differential_operators->D2DXDZ(f, outloc);
  }

  ASSERT2(result.getLocation() == outloc)

  return result;
}

FieldPerp Coordinates::Delp2(const FieldPerp& f, CELL_LOC outloc, bool useFFT) {
  TRACE("Coordinates::Delp2( FieldPerp )");

  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }

  ASSERT1(location == outloc)
  ASSERT1(f.getLocation() == outloc)

  if (localmesh->GlobalNx == 1 && localmesh->GlobalNz == 1) {
    // copy mesh, location, etc
    return f * 0;
  }
  ASSERT2(localmesh->xstart > 0) // Need at least one guard cell

  FieldPerp result{emptyFrom(f).setLocation(outloc)};

  int const jy = f.getIndex();
  result.setIndex(jy);

  if (useFFT) {
    int const ncz = localmesh->LocalNz;

    // Allocate memory
    auto ft = Matrix<dcomplex>(localmesh->LocalNx, ncz / 2 + 1);
    auto delft = Matrix<dcomplex>(localmesh->LocalNx, ncz / 2 + 1);

    // Take forward FFT
    for (int jx = 0; jx < localmesh->LocalNx; jx++) {
      rfft(&f(jx, 0), ncz, &ft(jx, 0));
    }

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
  }

  return result;
}

Field2D Coordinates::Laplace_par(const Field2D& f, CELL_LOC outloc) {
  ASSERT1(location == outloc || outloc == CELL_DEFAULT)

  return differential_operators->D2DY2(f, outloc) / g_22()
         + DDY(J() / g_22(), outloc) * DDY(f, outloc) / J();
}

Field3D Coordinates::Laplace_par(const Field3D& f, CELL_LOC outloc) {
  ASSERT1(location == outloc || outloc == CELL_DEFAULT)
  return differential_operators->D2DY2(f, outloc) / g_22()
         + differential_operators->DDY(J() / g_22(), outloc)
               * differential_operators->DDY(f, outloc) / J();
}

// Full Laplacian operator on scalar field

Field2D Coordinates::Laplace(const Field2D& f, CELL_LOC outloc,
                             const std::string& dfdy_boundary_conditions,
                             const std::string& dfdy_dy_region) {
  TRACE("Coordinates::Laplace( Field2D )");
  ASSERT1(location == outloc || outloc == CELL_DEFAULT)

  return G1() * DDX(f, outloc) + G2() * DDY(f, outloc)
         + g11() * differential_operators->D2DX2(f, outloc)
         + g22() * differential_operators->D2DY2(f, outloc)
         + 2.0 * g12()
               * differential_operators->D2DXDY(f, outloc, "DEFAULT", "RGN_NOBNDRY",
                                                dfdy_boundary_conditions, dfdy_dy_region);
}

Field3D Coordinates::Laplace(const Field3D& f, CELL_LOC outloc,
                             const std::string& dfdy_boundary_conditions,
                             const std::string& dfdy_dy_region) {
  TRACE("Coordinates::Laplace( Field3D )");
  ASSERT1(location == outloc || outloc == CELL_DEFAULT)

  return G1() * differential_operators->DDX(f, outloc)
         + G2() * differential_operators->DDY(f, outloc)
         + G3() * differential_operators->DDZ(f, outloc)
         + g11() * differential_operators->D2DX2(f, outloc)
         + g22() * differential_operators->D2DY2(f, outloc)
         + g33() * differential_operators->D2DZ2(f, outloc)
         + 2.0
               * (g12()
                      * differential_operators->D2DXDY(
                          f, outloc, "DEFAULT", "RGN_NOBNDRY", dfdy_boundary_conditions,
                          dfdy_dy_region)
                  + g13() * differential_operators->D2DXDZ(f, outloc)
                  + g23() * differential_operators->D2DYDZ(f, outloc));
}

// Full perpendicular Laplacian, in form of inverse of Laplacian operator in LaplaceXY
// solver
Field2D Coordinates::Laplace_perpXY(MAYBE_UNUSED(const Field2D& A),
                                    MAYBE_UNUSED(const Field2D& f)) const {
  TRACE("Coordinates::Laplace_perpXY( Field2D )");
#if not(BOUT_USE_METRIC_3D)
  Field2D result;
  result.allocate();
  for (auto i : result.getRegion(RGN_NOBNDRY)) {
    result[i] = 0.;

    // outer x boundary
    const auto outer_x_avg = [&i](const auto& f) { return 0.5 * (f[i] + f[i.xp()]); };
    const BoutReal outer_x_A = outer_x_avg(A);
    const BoutReal outer_x_J = outer_x_avg(J());
    const BoutReal outer_x_g11 = outer_x_avg(g11());
    const BoutReal outer_x_dx = outer_x_avg(dx());
    const BoutReal outer_x_value =
        outer_x_A * outer_x_J * outer_x_g11 / (J()[i] * outer_x_dx * dx()[i]);
    result[i] += outer_x_value * (f[i.xp()] - f[i]);

    // inner x boundary
    const auto inner_x_avg = [&i](const auto& f) { return 0.5 * (f[i] + f[i.xm()]); };
    const BoutReal inner_x_A = inner_x_avg(A);
    const BoutReal inner_x_J = inner_x_avg(J());
    const BoutReal inner_x_g11 = inner_x_avg(g11());
    const BoutReal inner_x_dx = inner_x_avg(dx());
    const BoutReal inner_x_value =
        inner_x_A * inner_x_J * inner_x_g11 / (J()[i] * inner_x_dx * dx()[i]);
    result[i] += inner_x_value * (f[i.xm()] - f[i]);

    // upper y boundary
    const auto upper_y_avg = [&i](const auto& f) { return 0.5 * (f[i] + f[i.yp()]); };
    const BoutReal upper_y_A = upper_y_avg(A);
    const BoutReal upper_y_J = upper_y_avg(J());
    const BoutReal upper_y_g_22 = upper_y_avg(g_22());
    const BoutReal upper_y_g23 = upper_y_avg(g23());
    const BoutReal upper_y_g_23 = upper_y_avg(g_23());
    const BoutReal upper_y_dy = upper_y_avg(dy());
    const BoutReal upper_y_value = -upper_y_A * upper_y_J * upper_y_g23 * upper_y_g_23
                                   / (upper_y_g_22 * J()[i] * upper_y_dy * dy()[i]);
    result[i] += upper_y_value * (f[i.yp()] - f[i]);

    // lower y boundary
    const auto lower_y_avg = [&i](const auto& f) { return 0.5 * (f[i] + f[i.ym()]); };
    const BoutReal lower_y_A = lower_y_avg(A);
    const BoutReal lower_y_J = lower_y_avg(J());
    const BoutReal lower_y_g_22 = lower_y_avg(g_22());
    const BoutReal lower_y_g23 = lower_y_avg(g23());
    const BoutReal lower_y_g_23 = lower_y_avg(g_23());
    const BoutReal lower_y_dy = lower_y_avg(dy());
    const BoutReal lower_y_value = -lower_y_A * lower_y_J * lower_y_g23 * lower_y_g_23
                                   / (lower_y_g_22 * J()[i] * lower_y_dy * dy()[i]);
    result[i] += lower_y_value * (f[i.ym()] - f[i]);
  }

  return result;
#else
  throw BoutException("Coordinates::Laplace_perpXY for 3D metric not implemented");
#endif
}

ChristoffelSymbols Coordinates::christoffel_symbols() const {
  if (christoffel_symbols_cache == nullptr) {
    auto ptr = std::make_unique<ChristoffelSymbols>(localmesh, differential_operators);
    ptr->CalculateChristoffelSymbols(*this);
    christoffel_symbols_cache = std::move(ptr);
  }
  return *christoffel_symbols_cache;
}

const Coordinates::FieldMetric& Coordinates::invSg() const {
  if (invSgCache == nullptr) {
    auto ptr = std::make_unique<FieldMetric>();
    (*ptr) = 1.0 / sqrt(g_22());
    invSgCache = std::move(ptr);
  }
  return *invSgCache;
}

const Coordinates::FieldMetric&
Coordinates::Grad2_par2_DDY_invSg(CELL_LOC outloc, const std::string& method) const {

  if (auto search = Grad2_par2_DDY_invSgCache.find(method);
      search != Grad2_par2_DDY_invSgCache.end()) {
    return *search->second;
  }
  invSg();

  // Communicate to get parallel slices
  localmesh->communicate(*invSgCache);
  invSgCache->applyParallelBoundary("parallel_neumann");

  // cache
  auto ptr = std::make_unique<FieldMetric>();
  *ptr = DDY(*invSgCache, outloc, method) * invSg();
  Grad2_par2_DDY_invSgCache[method] = std::move(ptr);
  return *Grad2_par2_DDY_invSgCache[method];
}

void Coordinates::checkCovariant() { covariantMetricTensor.check(localmesh->ystart); }

void Coordinates::checkContravariant() {
  contravariantMetricTensor.check(localmesh->ystart);
}

bool Coordinates::non_uniform() const { return non_uniform_; }
void Coordinates::setNon_uniform(bool non_uniform) { non_uniform_ = non_uniform; }

const FieldMetric& Coordinates::dx() const { return dx_; }
const FieldMetric& Coordinates::dy() const { return dy_; }
const FieldMetric& Coordinates::dz() const { return dz_; }

void Coordinates::setDx(FieldMetric dx) { dx_ = std::move(dx); }
void Coordinates::setDy(FieldMetric dy) { dy_ = std::move(dy); }
void Coordinates::setDz(FieldMetric dz) { dz_ = std::move(dz); }

void Coordinates::setDy(BoutReal value, int x, int y) { dy_(x, y) = value; }

void Coordinates::setD1_dx(FieldMetric d1_dx) { d1_dx_ = std::move(d1_dx); }
void Coordinates::setD1_dy(FieldMetric d1_dy) { d1_dy_ = std::move(d1_dy); }
void Coordinates::setD1_dz(FieldMetric d1_dz) { d1_dz_ = std::move(d1_dz); }

const FieldMetric& Coordinates::d1_dx() const { return d1_dx_; }
const FieldMetric& Coordinates::d1_dy() const { return d1_dy_; }
const FieldMetric& Coordinates::d1_dz() const { return d1_dz_; }

const FieldMetric& Coordinates::ShiftTorsion() const { return ShiftTorsion_; }

const FieldMetric& Coordinates::IntShiftTorsion() const { return IntShiftTorsion_; }

void Coordinates::setIntShiftTorsion(FieldMetric IntShiftTorsion) {
  IntShiftTorsion_ = std::move(IntShiftTorsion);
}

const MetricTensor::FieldMetric& Coordinates::g_11() const {
  return covariantMetricTensor.g11();
}
const MetricTensor::FieldMetric& Coordinates::g_22() const {
  return covariantMetricTensor.g22();
}
const MetricTensor::FieldMetric& Coordinates::g_33() const {
  return covariantMetricTensor.g33();
}
const MetricTensor::FieldMetric& Coordinates::g_12() const {
  return covariantMetricTensor.g12();
}
const MetricTensor::FieldMetric& Coordinates::g_13() const {
  return covariantMetricTensor.g13();
}
const MetricTensor::FieldMetric& Coordinates::g_23() const {
  return covariantMetricTensor.g23();
}

const MetricTensor::FieldMetric& Coordinates::g11() const {
  return contravariantMetricTensor.g11();
}
const MetricTensor::FieldMetric& Coordinates::g22() const {
  return contravariantMetricTensor.g22();
}
const MetricTensor::FieldMetric& Coordinates::g33() const {
  return contravariantMetricTensor.g33();
}
const MetricTensor::FieldMetric& Coordinates::g12() const {
  return contravariantMetricTensor.g12();
}
const MetricTensor::FieldMetric& Coordinates::g13() const {
  return contravariantMetricTensor.g13();
}
const MetricTensor::FieldMetric& Coordinates::g23() const {
  return contravariantMetricTensor.g23();
}

const FieldMetric& Coordinates::J() const { return J_; }

void Coordinates::setJ(FieldMetric J) {
  //TODO: Calculate J and check value is close
  J_ = std::move(J);
}

void Coordinates::setJ(BoutReal value, int x, int y) {
  //TODO: Calculate J and check value is close
  J_(x, y) = value;
}

const FieldMetric& Coordinates::Bxy() const { return Bxy_; }

void Coordinates::setBxy(FieldMetric Bxy) {
  //TODO: Calculate Bxy and check value is close
  Bxy_ = std::move(Bxy);
}

const FieldMetric& Coordinates::G1_11() const { return christoffel_symbols().G1_11(); }

const FieldMetric& Coordinates::G1_22() const { return christoffel_symbols().G1_22(); }

const FieldMetric& Coordinates::G1_33() const { return christoffel_symbols().G1_33(); }

const FieldMetric& Coordinates::G1_12() const { return christoffel_symbols().G1_12(); }

const FieldMetric& Coordinates::G1_13() const { return christoffel_symbols().G1_13(); }

const FieldMetric& Coordinates::G1_23() const { return christoffel_symbols().G1_23(); }

const FieldMetric& Coordinates::G2_11() const { return christoffel_symbols().G2_11(); }

const FieldMetric& Coordinates::G2_22() const { return christoffel_symbols().G2_22(); }

const FieldMetric& Coordinates::G2_33() const { return christoffel_symbols().G2_33(); }

const FieldMetric& Coordinates::G2_12() const { return christoffel_symbols().G2_12(); }

const FieldMetric& Coordinates::G2_13() const { return christoffel_symbols().G2_13(); }

const FieldMetric& Coordinates::G2_23() const { return christoffel_symbols().G2_23(); }

const FieldMetric& Coordinates::G3_11() const { return christoffel_symbols().G3_11(); }

const FieldMetric& Coordinates::G3_22() const { return christoffel_symbols().G3_22(); }

const FieldMetric& Coordinates::G3_33() const { return christoffel_symbols().G3_33(); }

const FieldMetric& Coordinates::G3_12() const { return christoffel_symbols().G3_12(); }

const FieldMetric& Coordinates::G3_13() const { return christoffel_symbols().G3_13(); }

const FieldMetric& Coordinates::G3_23() const { return christoffel_symbols().G3_23(); }

const FieldMetric& Coordinates::G1() const { return christoffel_symbols().G1(); }

const FieldMetric& Coordinates::G2() const { return christoffel_symbols().G2(); }

const FieldMetric& Coordinates::G3() const { return christoffel_symbols().G3(); }

void Coordinates::setG1(FieldMetric G1) { christoffel_symbols().setG1(G1); }
void Coordinates::setG2(FieldMetric G2) { christoffel_symbols().setG2(G2); }
void Coordinates::setG3(FieldMetric G3) { christoffel_symbols().setG3(G3); }

void Coordinates::applyToChristoffelSymbols(
    const std::function<const FieldMetric(const FieldMetric)>& function) {
  christoffel_symbols().map(function);
}

const MetricTensor& Coordinates::getContravariantMetricTensor() const {
  return contravariantMetricTensor;
}

const MetricTensor& Coordinates::getCovariantMetricTensor() const {
  return covariantMetricTensor;
}

void Coordinates::setContravariantMetricTensor(const MetricTensor& metric_tensor,
                                               const std::string& region,
                                               bool recalculate_staggered,
                                               bool force_interpolate_from_centre) {
  contravariantMetricTensor.setMetricTensor(metric_tensor);
  covariantMetricTensor.setMetricTensor(contravariantMetricTensor.inverse(region));
  recalculateAndReset(recalculate_staggered, force_interpolate_from_centre);
}

void Coordinates::setCovariantMetricTensor(const MetricTensor& metric_tensor,
                                           const std::string& region,
                                           bool recalculate_staggered,
                                           bool force_interpolate_from_centre) {
  covariantMetricTensor.setMetricTensor(metric_tensor);
  contravariantMetricTensor.setMetricTensor(covariantMetricTensor.inverse(region));
  recalculateAndReset(recalculate_staggered, force_interpolate_from_centre);
}

void Coordinates::applyToContravariantMetricTensor(
    const std::function<const FieldMetric(const FieldMetric)>& function) {
  contravariantMetricTensor.map(function);
}

void Coordinates::applyToCovariantMetricTensor(
    const std::function<const FieldMetric(const FieldMetric)>& function) {
  covariantMetricTensor.map(function);
}

void Coordinates::communicateChristoffelSymbolTerms() const {

  output_progress.write("\tCommunicating connection terms\n");

  auto tmp = G1_11(); // TODO: There must be a better way than this!
  communicate(tmp, G1_22(), G1_33(), G1_12(), G1_13(), G1_23(), G2_11(), G2_22(), G2_33(),
              G2_12(), G2_13(), G2_23(), G3_11(), G3_22(), G3_33(), G3_12(), G3_13(),
              G3_23(), G1(), G2(), G3());
}

void Coordinates::invalidateAndRecalculateCachedVariables() {
  zlength_cache.reset();
  Grad2_par2_DDY_invSgCache.clear();
  invSgCache.reset();
}