/**************************************************************************
 * Differential geometry
 * Calculates the covariant metric tensor, and christoffel symbol terms
 * given the contravariant metric tensor terms
 **************************************************************************/

#include "bout/christoffel_symbols.hxx"
#include "bout/coordinates_accessor.hxx"
#include "bout/field3d.hxx"
#include "bout/field_data.hxx"
#include "bout/g_values.hxx"
#include <bout/assert.hxx>
#include <bout/bout_types.hxx>
#include <bout/boutexception.hxx>
#include <bout/build_config.hxx>
#include <bout/build_defines.hxx>
#include <bout/constants.hxx>
#include <bout/coordinates.hxx>
#include <bout/dcomplex.hxx>
#include <bout/derivs.hxx>
#include <bout/fft.hxx>
#include <bout/field2d.hxx>
#include <bout/globals.hxx>
#include <bout/index_derivs_interface.hxx>
#include <bout/interpolation.hxx>
#include <bout/metric_tensor.hxx>
#include <bout/openmpwrap.hxx>
#include <bout/options.hxx>
#include <bout/output.hxx>
#include <bout/output_bout_types.hxx>
#include <bout/paralleltransform.hxx>
#include <bout/region.hxx>
#include <bout/sys/gettext.hxx>
#include <bout/sys/timer.hxx>
#include <bout/traits.hxx>
#include <bout/unused.hxx>
#include <bout/utils.hxx>
#include <bout/yboundary_regions.hxx>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <string>
#include <string_view>
#include <utility>

#include <fmt/format.h>

#include "parallel/fci.hxx"
#include "parallel/shiftedmetricinterp.hxx"

namespace {
// Extrapolate into boundaries (if requested) so that differential geometry
// terms can be interpolated if necessary
// Note: cannot use applyBoundary("free_o3") here because applyBoundary()
// would try to create a new Coordinates object since we have not finished
// initializing yet, leading to an infinite recursion.
// Also, here we interpolate for the boundary points at xstart/ystart and
// (xend+1)/(yend+1) instead of extrapolating.
template <class T, typename = bout::utils::EnableIfField<T>>
void fillGuards_impl(T& result, CELL_LOC location, const T& f, bool extrapolate_x,
                     bool extrapolate_y, bool no_extra_interpolate = false) {
  const auto* localmesh = result.getMesh();

  const auto zstart = bout::build::use_metric_3d ? localmesh->zstart : 0;
  const auto zend = bout::build::use_metric_3d ? localmesh->zend : 0;

  for (auto& newbndry : localmesh->getBoundaries()) {
    auto* bndry = newbndry->getLegacyPointer();
    if ((extrapolate_x and bndry->bx != 0) or (extrapolate_y and bndry->by != 0)) {
      // Can use no_extra_interpolate argument to skip the extra interpolation when we
      // want to extrapolate the Christoffel symbol terms which come from derivatives so
      // don't have the extra point set already
      const bool interpolating_into_x_boundary =
          ((location == CELL_XLOW) and (bndry->bx > 0));
      const bool interpolating_into_y_boundary =
          ((location == CELL_YLOW) and (bndry->by > 0));
      const int extrap_start =
          (not no_extra_interpolate)
                  and (interpolating_into_x_boundary or interpolating_into_y_boundary)
              ? 1
              : 0;
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
          for (int zi = zstart; zi <= zend; ++zi) {
            result(bndry->x, bndry->y, zi) =
                (9.
                     * (f(bndry->x - bndry->bx, bndry->y - bndry->by, zi)
                        + f(bndry->x, bndry->y, zi))
                 - f(bndry->x - (2 * bndry->bx), bndry->y - (2 * bndry->by), zi)
                 - f(bndry->x + bndry->bx, bndry->y + bndry->by, zi))
                / 16.;
          }
        }

        // set boundary guard cells
        if ((bndry->bx != 0 && localmesh->GlobalNx - (2 * bndry->width) >= 3)
            || (bndry->by != 0
                && localmesh->GlobalNy - (localmesh->numberOfYBoundaries() * bndry->width)
                       >= 3)) {
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
            const int xi = bndry->x + (i * bndry->bx);
            const int yi = bndry->y + (i * bndry->by);
            for (int zi = zstart; zi <= zend; ++zi) {
              result(xi, yi, zi) =
                  (3.0 * result(xi - bndry->bx, yi - bndry->by, zi))
                  - (3.0 * result(xi - (2 * bndry->bx), yi - (2 * bndry->by), zi))
                  + result(xi - (3 * bndry->bx), yi - (3 * bndry->by), zi);
            }
          }
        } else {
          // not enough grid points to extrapolate, set equal to last grid point
          for (int i = extrap_start; i < bndry->width; i++) {
            for (int zi = zstart; zi <= zend; ++zi) {
              result(bndry->x + (i * bndry->bx), bndry->y + (i * bndry->by), zi) =
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
        for (int k = 0; k < f.getNz(); ++k) {
          result(i, j, k) = BoutNaN;
          result(i, localmesh->LocalNy - 1 - j, k) = BoutNaN;
          result(localmesh->LocalNx - 1 - i, j, k) = BoutNaN;
          result(localmesh->LocalNx - 1 - i, localmesh->LocalNy - 1 - j, k) = BoutNaN;
        }
      }
    }
  }
#endif
}

/// Interpolate a Field to a new CELL_LOC with interp_to.
/// Communicates to set internal guard cells.
/// Boundary guard cells are set by extrapolating from the grid, like
/// 'free_o3' boundary conditions
/// Corner guard cells are set to BoutNaN
template <class T, typename = bout::utils::EnableIfField<T>>
auto interpolateAndExtrapolate(const T& f_, CELL_LOC location, bool extrapolate_x,
                               bool extrapolate_y, bool no_extra_interpolate,
                               ParallelTransform* pt_ = nullptr) -> T {

  Mesh* localmesh = f_.getMesh();
  T f = f_;

  // If input f is member of the Coordinates we are currently constructing, it will not
  // have Coordinates and needs to use the passed-in ParallelTransform.
  // Otherwise, if input f is from Coordinates at a different location, it will have its own
  // Coordinates, and we should use its ParallelTransform
  ParallelTransform* pt_f =
      (not bout::build::use_metric_3d) or f.getCoordinates() == nullptr
          ? pt_
          : &f.getCoordinates()->getParallelTransform();

  if (f.getDirectionY() != YDirectionType::Standard) {
    if (pt_f->canToFromFieldAligned()) {
      f = pt_f->fromFieldAligned(f);
    } else {
      f.setDirectionY(YDirectionType::Standard);
    }
  }

  T result;
  if (location == CELL_YLOW and f.getLocation() != CELL_YLOW) {
    auto f_aligned = pt_f->toFieldAligned(f, "RGN_NOX");
    result = interp_to(f_aligned, location, "RGN_NOBNDRY");
    const auto* result_coords = result.getCoordinates();
    ParallelTransform* pt_result =
        result_coords == nullptr ? pt_ : &result_coords->getParallelTransform();

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
  localmesh->communicate_no_slices(result);

  fillGuards_impl(result, location, f, extrapolate_x, extrapolate_y,
                  no_extra_interpolate);

  return result;
}

// If the CELL_CENTRE variable was read, the staggered version is required to
// also exist for consistency
void checkStaggeredGet(Mesh* mesh, const std::string& name, const std::string& suffix) {
  if (mesh->sourceHasVar(name) != mesh->sourceHasVar(name + suffix)) {
    throw BoutException("Attempting to read staggered fields from grid, but '{}'"
                        " is not present in both CELL_CENTRE and staggered versions.",
                        name);
  }
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

Coordinates::Coordinates(Mesh* mesh, FieldMetric dx, FieldMetric dy, FieldMetric dz,
                         FieldMetric J, FieldMetric Bxy, FieldMetric g11, FieldMetric g22,
                         FieldMetric g33, FieldMetric g12, FieldMetric g13,
                         FieldMetric g23, FieldMetric g_11, FieldMetric g_22,
                         FieldMetric g_33, FieldMetric g_12, FieldMetric g_13,
                         FieldMetric g_23, FieldMetric ShiftTorsion,
                         FieldMetric IntShiftTorsion)
    : nz(mesh->LocalNz), localmesh(mesh), location(CELL_CENTRE), dx_(std::move(dx)),
      dy_(std::move(dy)), dz_(std::move(dz)), ShiftTorsion_(std::move(ShiftTorsion)),
      IntShiftTorsion_(std::move(IntShiftTorsion)),
      transform(bout::utils::make_unique<ParallelTransformIdentity>(*localmesh)),
      contravariantMetricTensor(std::move(g11), std::move(g22), std::move(g33),
                                std::move(g12), std::move(g13), std::move(g23)),
      covariantMetricTensor(std::move(g_11), std::move(g_22), std::move(g_33),
                            std::move(g_12), std::move(g_13), std::move(g_23)),
      jacobian_cache(std::make_unique<FieldMetric>(std::move(J))), Bxy_(std::move(Bxy)) {}

Coordinates::Coordinates(Mesh* mesh, Options* options)
    : nz(mesh->LocalNz), localmesh(mesh), localoptions(options), location(CELL_CENTRE),
      dx_(1., mesh), dy_(1., mesh), dz_(1., mesh), d1_dx_(mesh), d1_dy_(mesh),
      d1_dz_(mesh), ShiftTorsion_(0.0, mesh), IntShiftTorsion_(0.0, mesh),
      contravariantMetricTensor(1., 1., 1., 0, 0, 0, mesh),
      covariantMetricTensor(1., 1., 1., 0, 0, 0, mesh), Bxy_(1., mesh) {

  readFromMesh(options, "");

  if (transform != nullptr and not transform->canToFromFieldAligned()) {
    // Read parallel metrics from gridfile for FCI
    readParallelMetricComponents();
  }
}

Coordinates::Coordinates(Mesh* mesh, Options* options, const CELL_LOC loc,
                         const Coordinates* coords_in, bool force_interpolate_from_centre)
    : nz(mesh->LocalNz), localmesh(mesh), localoptions(options), location(loc),
      dx_(1., mesh), dy_(1., mesh), dz_(1., mesh), d1_dx_(mesh), d1_dy_(mesh),
      d1_dz_(mesh), ShiftTorsion_(0.0, mesh), IntShiftTorsion_(0.0, mesh),
      contravariantMetricTensor(1., 1., 1., 0, 0, 0, mesh),
      covariantMetricTensor(1., 1., 1., 0, 0, 0, mesh), Bxy_(1., mesh) {

  const std::string suffix = getLocationSuffix(location);

  if (not force_interpolate_from_centre and mesh->sourceHasVar("dx" + suffix)) {
    readFromMesh(options, suffix);
  } else {
    // Interpolate fields from coords_in

    if (isUniform(coords_in->dz_)) {
      dz_ = coords_in->dz_;
      dz_.setLocation(location);
    } else {
      throw BoutException(
          "We are asked to transform dz to get dz before we have a transform, which "
          "might require dz!\nPlease provide a dz for the staggered quantity!");
    }
    setParallelTransform(options);

    auto interpField = [&, this](const FieldMetric& f) -> FieldMetric {
      return interpolateAndExtrapolate(f, location, true, true, false, transform.get());
    };

    dx_ = interpField(coords_in->dx_);
    dy_ = interpField(coords_in->dy_);
    // not really needed - we have used dz already ...
    dz_ = interpField(coords_in->dz_);

    setMetricTensor(coords_in->getContravariantMetricTensor(),
                    coords_in->getCovariantMetricTensor());

    contravariantMetricTensor.map(interpField);
    covariantMetricTensor.map(interpField);

    // Check input metrics
    checkContravariant();
    checkCovariant();

    setJ(interpField(coords_in->J()));
    setBxy(interpField(coords_in->Bxy()));

    bout::checkFinite(J(), "The Jacobian", "RGN_NOCORNERS");
    bout::checkPositive(J(), "The Jacobian", "RGN_NOCORNERS");
    bout::checkFinite(Bxy(), "Bxy", "RGN_NOCORNERS");
    bout::checkPositive(Bxy(), "Bxy", "RGN_NOCORNERS");

    ShiftTorsion_ = interpField(coords_in->ShiftTorsion_);

    if (mesh->IncIntShear) {
      IntShiftTorsion_ = interpField(coords_in->IntShiftTorsion_);
    }
    if (not transform->canToFromFieldAligned()) {
      // Read parallel metrics from gridfile for FCI
      readParallelMetricComponents();
    }
  }
}

void Coordinates::readFromMesh(Options* options, const std::string& suffix) {
  if (options == nullptr) {
    options = Options::getRoot()->getSection("mesh");
  }

  // Note: If boundary cells were not loaded from the grid file, use
  // 'interpolateAndExtrapolate' to set them. Ensures that derivatives are
  // smooth at all the boundaries.

  const bool extrapolate_x =
      (*options)["extrapolate_x"].withDefault(not localmesh->sourceHasXBoundaryGuards());
  const bool extrapolate_y =
      (*options)["extrapolate_y"].withDefault(not localmesh->sourceHasYBoundaryGuards());

  if (extrapolate_x) {
    output_warn.write(_("WARNING: extrapolating input mesh quantities into x-boundary "
                        "cells. Set option extrapolate_x=false to disable this.\n"));
  }

  if (extrapolate_y) {
    output_warn.write(_("WARNING: extrapolating input mesh quantities into y-boundary "
                        "cells. Set option extrapolate_y=false to disable this.\n"));
  }

  // Some helper functions that let us avoid passing the same arguments every time

  // Ensure that we have an unaligned field, transforming if necessary
  auto ensuredUnaligned = [this](FieldMetric& field) -> FieldMetric {
    if (field.getDirectionY() == YDirectionType::Aligned
        and transform->canToFromFieldAligned()) {
      return transform->fromFieldAligned(field);
    }
    return field.setDirectionY(YDirectionType::Standard);
  };

  // Read from the mesh and transform from field aligned if required
  auto readField = [this, &suffix, ensuredUnaligned](const std::string& name,
                                                     BoutReal def_value) -> FieldMetric {
    checkStaggeredGet(localmesh, name, suffix);

    FieldMetric field{localmesh};
    localmesh->get(field, name + suffix, def_value, false, location);
    return ensuredUnaligned(field);
  };

  // Just fill in the guard cells (via extrapolation) for an existing field
  auto fillGuards = [&, this](FieldMetric& field) {
    // Passing `field` in here twice is gross but ok because the second argument
    // is only used when interpolating, and we're not interpolating here
    fillGuards_impl(field, location, field, extrapolate_x, extrapolate_y);
    return field;
  };

  // Read the field, transform if required, and fill in the guards
  auto readAndFillGuards = [&](const std::string& name,
                               BoutReal def_value) -> FieldMetric {
    auto field = readField(name, def_value);
    fillGuards(field);
    return field;
  };

  {
    // dz
    auto& options = Options::root();
    const bool has_zperiod = options.isSet("zperiod");
    const auto zmin = has_zperiod ? 0.0 : options["ZMIN"].withDefault(0.0);
    const auto zmax = has_zperiod ? 1.0 / options["zperiod"].withDefault(1.0)
                                  : options["ZMAX"].withDefault(1.0);

    const auto default_dz = (zmax - zmin) * TWOPI / nz;

    // We can't use the helper functions here because in 3D we might (always?)
    // need to transform from field aligned -- which requires dz! So we have to
    // read it "plain" first...
    localmesh->get(dz_, "dz" + suffix, default_dz, false, location);
  }

  // ...then we can set the transform (required early for differentiation)...
  setParallelTransform(options);

  // ...and finally we can transform/interpolate/fill in guards
  setDz(interpolateAndExtrapolate(dz_, location, extrapolate_x, extrapolate_y, false,
                                  transform.get()));

  // everything else from this point on can use our helper functions
  setDx(readAndFillGuards("dx", 1.0), localmesh->periodicX);

  setDy(readAndFillGuards("dy", 1.0));

  // Diagonal components of metric tensor g^{ij} (default to 1)
  contravariantMetricTensor = ContravariantMetricTensor{
      readAndFillGuards("g11", 1.0), readAndFillGuards("g22", 1.0),
      readAndFillGuards("g33", 1.0), readAndFillGuards("g12", 0.0),
      readAndFillGuards("g13", 0.0), readAndFillGuards("g23", 0.0),
  };

  // Check input metrics
  checkContravariant();

  // Find covariant metric components
  auto covariant_component_names = {"g_11", "g_22", "g_33", "g_12", "g_13", "g_23"};
  auto source_has_component = [&suffix, this](const std::string& name) {
    return localmesh->sourceHasVar(name + suffix);
  };

  // Check if any of the components are present, all of them are present
  if (const bool all_present =
          std::ranges::all_of(covariant_component_names, source_has_component);
      std::ranges::any_of(covariant_component_names, source_has_component)
      and all_present) {
    covariantMetricTensor = CovariantMetricTensor{
        readAndFillGuards("g_11", 1.0), readAndFillGuards("g_22", 1.0),
        readAndFillGuards("g_33", 1.0), readAndFillGuards("g_12", 0.0),
        readAndFillGuards("g_13", 0.0), readAndFillGuards("g_23", 0.0),
    };

    output_warn.write("\tWARNING! Covariant components of metric tensor set manually. "
                      "Contravariant components NOT recalculated\n");
  } else {
    if (not all_present) {
      output_warn.write("Not all covariant components of metric tensor found. "
                        "Calculating all from the contravariant tensor\n");
    }
    covariantMetricTensor = contravariantMetricTensor.inverse("RGN_NOCORNERS", false);

    // More robust to extrapolate derived quantities directly, rather than
    // deriving from extrapolated covariant metric components
    covariantMetricTensor.map(fillGuards);
  }

  // Check covariant metrics
  checkCovariant();

  // Attempt to read J from the grid file
  auto Jcalc = J();
  FieldMetric J_temp{localmesh};
  if (localmesh->get(J_temp, "J" + suffix, 0.0, false) != 0) {
    output_warn.write(
        "\tWARNING: Jacobian 'J' not found. Calculating from metric tensor\n");
    localmesh->communicate_no_slices(*jacobian_cache);
  } else {
    checkStaggeredGet(localmesh, "J", suffix);
    *jacobian_cache = ensuredUnaligned(J_temp);
    fillGuards(*jacobian_cache);

    // Compare calculated and loaded values
    output_warn.write("\tMaximum difference in J is {:e}\n", max(abs(J() - Jcalc)));

    localmesh->communicate_no_slices(*jacobian_cache);
  }

  // Check jacobian
  bout::checkFinite(J(), "J" + suffix, "RGN_NOCORNERS");
  bout::checkPositive(J(), "J" + suffix, "RGN_NOCORNERS");
  if (min(abs(J())) < 1.0e-10) {
    throw BoutException("\tERROR: Jacobian{:s} becomes very small\n", suffix);
  }

  // Attempt to read Bxy from the grid file
  const FieldMetric Bcalc = sqrt(g_22()) / J();
  if (localmesh->get(Bxy_, "Bxy" + suffix, 0.0, false) != 0) {
    output_warn.write("\tWARNING: Magnitude of B field 'Bxy' not found. Calculating from "
                      "metric tensor\n");
    Bxy_ = interpolateAndExtrapolate(Bcalc, location, extrapolate_x, extrapolate_y, false,
                                     transform.get());
    localmesh->communicate_no_slices(Bxy_);
  } else {
    checkStaggeredGet(localmesh, "Bxy", suffix);
    Bxy_ = ensuredUnaligned(Bxy_);
    fillGuards(Bxy_);
    output_warn.write("\tMaximum difference in Bxy is {:e}\n", max(abs(Bxy_ - Bcalc)));
  }

  // Check Bxy
  bout::checkFinite(Bxy(), "Bxy" + suffix, "RGN_NOCORNERS");
  bout::checkPositive(Bxy(), "Bxy" + suffix, "RGN_NOCORNERS");

  if (not localmesh->sourceHasVar("ShiftTorsion" + suffix)) {
    output_warn.write("\tWARNING: No Torsion specified for zShift. "
                      "Derivatives may not be correct\n");
  } else {
    ShiftTorsion_ = readAndFillGuards("ShiftTorsion", 0.0);
  }

  if (localmesh->IncIntShear) {
    if (not localmesh->sourceHasVar("IntShiftTorsion" + suffix)) {
      output_warn.write("\tWARNING: No Integrated torsion specified\n");
    } else {
      IntShiftTorsion_ = readAndFillGuards("IntShiftTorsion", 0.0);
    }
  }
}

namespace bout {
// Get a unique name for a field based on the sign/magnitude of the offset
std::string parallelSliceFieldName(std::string_view field, int offset) {
  using namespace std::string_view_literals;
  const std::string_view direction = (offset > 0) ? "forward"sv : "backward"sv;
  // We only have a suffix for parallel slices beyond the first
  // This is for backwards compatibility
  const std::string slice_suffix =
      (std::abs(offset) > 1) ? fmt::format("_{}", std::abs(offset)) : "";
  return fmt::format("{}_{}{}", direction, field, slice_suffix);
};
} // namespace bout

namespace {
void load_parallel_metric_component(std::string_view name, bout::FieldMetric& component,
                                    int offset) {
  Mesh* mesh = component.getMesh();
  bout::FieldMetric tmp{mesh};
  const auto pname = bout::parallelSliceFieldName(name, offset);
  if (mesh->get(tmp, pname, 0.0, false) != 0) {
    throw BoutException("Could not read {:s} from grid file!\n"
                        "  Fix it up with `zoidberg-update-parallel-metrics <grid>`",
                        pname);
  }
  if (!component.hasParallelSlices()) {
    component.splitParallelSlices();
    component.disallowCalcParallelSlices();
    component.resetRegionParallel(true);
  }
  auto& pcom = component.ynext(offset);
  pcom.allocate();
  BOUT_FOR(i, component.getRegion("RGN_NOBNDRY")) { pcom[i.yp(offset)] = tmp[i]; }
}
} // namespace

void Coordinates::readParallelMetricComponents() {
  if (not bout::build::use_metric_3d) {
    return;
  }
  output_info.write("\tLoading parallel metrics\n");
  const FieldMetric JB0 = J() * Bxy();
  jacobian_cache->splitParallelSlices();
  jacobian_cache->disallowCalcParallelSlices();
  jacobian_cache->resetRegionParallel(true);
  for (int i = 1; i <= localmesh->ystart; ++i) {
    auto read_offset = [i](std::string_view name, FieldMetric& component) {
      load_parallel_metric_component(name, component, -i);
      load_parallel_metric_component(name, component, i);
    };

    read_offset("g_11", covariantMetricTensor.g11_m);
    read_offset("g_22", covariantMetricTensor.g22_m);
    read_offset("g_33", covariantMetricTensor.g33_m);
    read_offset("g_13", covariantMetricTensor.g13_m);
    read_offset("g11", contravariantMetricTensor.g11_m);
    read_offset("g22", contravariantMetricTensor.g22_m);
    read_offset("g33", contravariantMetricTensor.g33_m);
    read_offset("g13", contravariantMetricTensor.g13_m);
    read_offset("dy", dy_);
    read_offset("Bxy", Bxy_);

    jacobian_cache->ynext(i).allocate();
    jacobian_cache->ynext(-i).allocate();
    BOUT_FOR(j, JB0.getRegion("RGN_NOBNDRY")) {
      jacobian_cache->ynext(i)[j.yp(i)] = JB0[j] / Bxy_.ynext(i)[j.yp(i)];
      jacobian_cache->ynext(-i)[j.yp(-i)] = JB0[j] / Bxy_.ynext(-i)[j.yp(-i)];
    }
  }
}

void Coordinates::outputVars(Options& output_options) {
  const Timer time("io");
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

  if (g_values_cache != nullptr) {
    // If we haven't used them yet, then presumably we don't actually need them,
    // so let's not compute them now.
    // Also, previously, we were explicitly setting the G-values to NaN for FCI
    // instead of computing them (presumably because these are both difficult to
    // compute accurately for FCI and unneeded in FCI models in practice).
    output_options["G1" + loc_string].force(G1(), "Coordinates");
    output_options["G2" + loc_string].force(G2(), "Coordinates");
    output_options["G3" + loc_string].force(G3(), "Coordinates");
  }

  getParallelTransform().outputVars(output_options);
}

const Field2D& Coordinates::zlength() const {
  BOUT_OMP(critical)
  if (not zlength_cache) {
    zlength_cache = std::make_unique<Field2D>(0., localmesh);

#if BOUT_USE_METRIC_3D
    BOUT_FOR_SERIAL(i, dz().getRegion("RGN_ALL")) { (*zlength_cache)[i] += dz()[i]; }
#else
    (*zlength_cache) = dz_ * nz;
#endif
  }

  return *zlength_cache;
}

void Coordinates::setDx(FieldMetric dx, const bool communicate) {
  if (min(abs(dx)) < 1e-8) {
    throw BoutException("dx magnitude less than 1e-8");
  }
  dx_ = std::move(dx);
  invalidateCellGeometryCaches();
  invalidateAccessorCache();
  if (communicate) {
    localmesh->communicate_no_slices(dx_);
  }
}

void Coordinates::setDy(FieldMetric dy, const bool communicate) {
  if (min(abs(dy)) < 1e-8) {
    throw BoutException("dy magnitude less than 1e-8");
  }
  dy_ = std::move(dy);
  invalidateCellGeometryCaches();
  invalidateAccessorCache();
  if (communicate) {
    localmesh->communicate_no_slices(dy_);
  }
}

void Coordinates::setDz(FieldMetric dz, const bool communicate) {
  if (min(abs(dz)) < 1e-8) {
    throw BoutException("dz magnitude less than 1e-8");
  }
  dz_ = std::move(dz);
  zlength_cache.reset();
  invalidateCellGeometryCaches();
  invalidateAccessorCache();
  if (communicate) {
    localmesh->communicate_no_slices(dz_);
  }
}

void Coordinates::recalculateAndReset(bool recalculate_staggered,
                                      bool force_interpolate_from_centre) {

  // Check input metrics
  checkContravariant();
  checkCovariant();

  christoffel_symbols_cache.reset();
  g_values_cache.reset();

  correctionForNonUniformMeshes(force_interpolate_from_centre);

  if (location == CELL_CENTRE && recalculate_staggered) {
    // Re-calculate interpolated Coordinates at staggered locations
    localmesh->recalculateStaggeredCoordinates();
  }

  zlength_cache.reset();
  Grad2_par2_DDY_invSgCache.clear();
  invSgCache.reset();
  invalidateCellGeometryCaches();
  invalidateAccessorCache();
}

void Coordinates::correctionForNonUniformMeshes(bool force_interpolate_from_centre) {
  OPTION(Options::getRoot(), non_uniform_, true);

  FieldMetric d2x(localmesh);
  FieldMetric d2y(localmesh);

  // Read correction for non-uniform meshes
  const std::string suffix = getLocationSuffix(location);

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
    d1_dx_ = bout::derivatives::index::DDX(FieldMetric{1. / dx()}); // d/di(1/dx)

    localmesh->communicate_no_slices(d1_dx_);
    d1_dx_ =
        interpolateAndExtrapolate(d1_dx_, location, true, true, true, transform.get());
  } else {
    d2x.setLocation(location);
    // set boundary cells if necessary
    d2x = interpolateAndExtrapolate(d2x, location, extrapolate_x, extrapolate_y, false,
                                    transform.get());

    d1_dx_ = -d2x / (dx() * dx());
  }

  if (localmesh->get(d2y, "d2y" + suffix, 0.0, false, location) != 0) {
    output_warn.write("\tWARNING: differencing quantity 'd2y' not found. "
                      "Calculating from dy\n");
    d1_dy_ = DDY(1. / dy()); // d/di(1/dy)

    localmesh->communicate_no_slices(d1_dy_);
    d1_dy_ =
        interpolateAndExtrapolate(d1_dy_, location, true, true, true, transform.get());
  } else {
    d2y.setLocation(location);
    // set boundary cells if necessary
    d2y = interpolateAndExtrapolate(d2y, location, extrapolate_x, extrapolate_y, false,
                                    transform.get());

    d1_dy_ = -d2y / (dy() * dy());
  }

  if (bout::build::use_metric_3d) {
    FieldMetric d2z(localmesh); // d^2 x / d i^2
    if (localmesh->get(d2z, "d2z" + suffix, 0.0, false, location) != 0) {
      output_warn.write("\tWARNING: differencing quantity 'd2z' not found. "
                        "Calculating from dz\n");
      d1_dz_ = bout::derivatives::index::DDZ(FieldMetric{1. / dz()});
      localmesh->communicate_no_slices(d1_dz_);
      d1_dz_ =
          interpolateAndExtrapolate(d1_dz_, location, true, true, true, transform.get());
    } else {
      d2z.setLocation(location);
      // set boundary cells if necessary
      d2z = interpolateAndExtrapolate(d2z, location, extrapolate_x, extrapolate_y, false,
                                      transform.get());

      d1_dz_ = -d2z / (dz() * dz());
    }
  } else {
    d1_dz_ = 0;
  }

  localmesh->communicate_no_slices(d1_dx_, d1_dy_, d1_dz_);
}

Coordinates::FieldMetric Coordinates::recalculateJacobian() const {
  // calculate Jacobian using g^-1 = det[g^ij], J = sqrt(g)
  const FieldMetric g_matrix = g11() * g22() * g33() + 2.0 * g12() * g13() * g23()
                               - g11() * g23() * g23() - g22() * g13() * g13()
                               - g33() * g12() * g12();

  bout::checkPositive(g_matrix, "The determinant of g^ij", "RGN_NOBNDRY");

  return 1. / sqrt(g_matrix);
}

Coordinates::FieldMetric Coordinates::recalculateBxy() const {
  return sqrt(g_22()) / J();
}

namespace {
// Utility function for fixing up guard cells of zShift
void fixZShiftGuards(Field2D& zShift) {
  auto* localmesh = zShift.getMesh();

  // extrapolate into boundary guard cells if necessary
  zShift = interpolateAndExtrapolate(zShift, zShift.getLocation(),
                                     not localmesh->sourceHasXBoundaryGuards(),
                                     not localmesh->sourceHasYBoundaryGuards(), false);

  // make sure zShift has been communicated
  localmesh->communicate_no_slices(zShift);

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

  auto ptoptions = options->getSection("paralleltransform");

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
    std::string suffix = getLocationSuffix(location);
    if (localmesh->sourceHasVar("dx" + suffix)) {
      // Grid file has variables at this location, so should be able to read
      checkStaggeredGet(localmesh, "zShift", suffix);
      if (localmesh->get(zShift, "zShift" + suffix, 0.0, false, location)) {
        // No zShift variable. Try qinty in BOUT grid files
        if (localmesh->get(zShift, "qinty" + suffix, 0.0, false, location)) {
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
      if (localmesh->get(zShift_centre, "zShift", 0.0, false)) {
        // No zShift variable. Try qinty in BOUT grid files
        if (localmesh->get(zShift_centre, "qinty", 0.0, false)) {
          // Failed to find either variable, cannot use ShiftedMetric
          throw BoutException("Could not read zShift from grid file");
        }
      }

      fixZShiftGuards(zShift_centre);

      zShift = interpolateAndExtrapolate(zShift_centre, location, true, true, false,
                                         transform.get());
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
    throw BoutException(_f("Unrecognised paralleltransform option '{}'.\n"
                           "Valid choices are 'identity', 'shifted', 'fci'"),
                        ptstr);
  }
}

const ChristoffelSymbols& Coordinates::christoffel_symbols() const {
  BOUT_OMP_SAFE(critical(christoffel_symbols_cache))
  {
    if (christoffel_symbols_cache == nullptr) {
      christoffel_symbols_cache = std::make_unique<ChristoffelSymbols>(*this);
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

      christoffel_symbols_cache->map([this](const FieldMetric& component) {
        return interpolateAndExtrapolate(component, location, true, true, false,
                                         transform.get());
      });
    }
  }
  return *christoffel_symbols_cache;
}

GValues& Coordinates::g_values() const {
  BOUT_OMP_SAFE(critical(g_values_cache))
  {
    if (g_values_cache == nullptr) {
      g_values_cache = std::make_unique<GValues>(*this);
      g_values_cache->map([this](const FieldMetric& component) {
        return interpolateAndExtrapolate(component, location, true, true, true,
                                         transform.get());
      });
    }
  }
  return *g_values_cache;
}

const Coordinates::FieldMetric& Coordinates::invSg() const {
  BOUT_OMP_SAFE(critical(invSg_cache))
  {
    if (invSgCache == nullptr) {
      auto ptr = std::make_unique<Coordinates::FieldMetric>();
      (*ptr) = 1.0 / sqrt(g_22());
      invSgCache = std::move(ptr);
    }
  }
  return *invSgCache;
}

const Coordinates::FieldMetric&
Coordinates::Grad2_par2_DDY_invSg(CELL_LOC outloc, const std::string& method) const {
  const FieldMetric* result{nullptr};
  BOUT_OMP_SAFE(critical(Grad2_par2_DDY_invSg_cache))
  {
    if (auto search = Grad2_par2_DDY_invSgCache.find(method);
        search != Grad2_par2_DDY_invSgCache.end()) {
      result = search->second.get();
    } else {
      if (invSgCache == nullptr) {
        auto ptr = std::make_unique<Coordinates::FieldMetric>();
        (*ptr) = 1.0 / sqrt(g_22());
        invSgCache = std::move(ptr);
      }

      // Communicate to get parallel slices
      localmesh->communicate(*invSgCache);
      invSgCache->applyParallelBoundary("parallel_neumann_o2");

      auto ptr = std::make_unique<Coordinates::FieldMetric>();
      *ptr = DDY(*invSgCache, outloc, method) * (*invSgCache);
      result = ptr.get();
      Grad2_par2_DDY_invSgCache[method] = std::move(ptr);
    }
  }
  return *result;
}

void Coordinates::checkCovariant() { covariantMetricTensor.check(localmesh->ystart); }

void Coordinates::checkContravariant() {
  contravariantMetricTensor.check(localmesh->ystart);
}

void Coordinates::invalidateCellGeometryCaches() {
  _g_22_ylow.reset();
  _g_22_yhigh.reset();
  _cell_area_xlow.reset();
  _cell_area_xhigh.reset();
  _cell_area_ylow.reset();
  _cell_area_yhigh.reset();
  _cell_area_zlow.reset();
  _cell_area_zhigh.reset();
  _cell_volume.reset();
}

void Coordinates::invalidateAccessorCache() const { CoordinatesAccessor::clear(this); }

void Coordinates::invalidateJacobianCaches() {
  g_values_cache.reset();
  invalidateCellGeometryCaches();
  invalidateAccessorCache();
}

void Coordinates::invalidateMetricCaches() {
  christoffel_symbols_cache.reset();
  g_values_cache.reset();
  Grad2_par2_DDY_invSgCache.clear();
  invSgCache.reset();
  jacobian_cache.reset();
  invalidateCellGeometryCaches();
  invalidateAccessorCache();
}

const Coordinates::FieldMetric& Coordinates::J() const {
  BOUT_OMP_SAFE(critical(jacobian_cache))
  {
    if (jacobian_cache == nullptr) {
      jacobian_cache = std::make_unique<FieldMetric>(recalculateJacobian());
    }
  }
  return *jacobian_cache;
}

void Coordinates::setJ(const FieldMetric& J, const bool communicate) {
  bout::checkFinite(J, "J", "RGN_NOCORNERS");
  bout::checkPositive(J, "J", "RGN_NOCORNERS");

  //TODO: Calculate J and check value is close
  invalidateJacobianCaches();
  jacobian_cache = std::make_unique<FieldMetric>(J);
  if (communicate) {
    localmesh->communicate_no_slices(*jacobian_cache);
  }
}

const Coordinates::FieldMetric& Coordinates::g_22_ylow() const {
  if (_g_22_ylow.has_value()) {
    return *_g_22_ylow;
  }
  BOUT_OMP_SAFE(critical)
  {
    if (!_g_22_ylow.has_value()) {
      _g_22_ylow.emplace(emptyFrom(g_22()));
      if (Bxy().isFci()) {
        if (localmesh->get(_g_22_ylow.value(), "g_22_cell_ylow", 0.0, false) != 0) {
          throw BoutException("The grid file does not contain `g_22_cell_ylow`.");
        }
      } else {
        ASSERT0(localmesh->ystart > 0);
        BOUT_FOR(i, g_22().getRegion("RGN_NOY")) {
          _g_22_ylow.value()[i] =
              SQ(0.5 * (std::sqrt(g_22()[i]) + std::sqrt(g_22()[i.ym()])));
        }
      }
    }
  }
  return g_22_ylow();
}

const Coordinates::FieldMetric& Coordinates::g_22_yhigh() const {
  if (_g_22_yhigh.has_value()) {
    return *_g_22_yhigh;
  }
  BOUT_OMP_SAFE(critical)
  {
    if (!_g_22_yhigh.has_value()) {
      _g_22_yhigh.emplace(emptyFrom(g_22()));
      if (Bxy().isFci()) {
        if (localmesh->get(_g_22_yhigh.value(), "g_22_cell_yhigh", 0.0, false) != 0) {
          throw BoutException("The grid file does not contain `g_22_cell_yhigh`.");
        }
      } else {
        ASSERT0(localmesh->ystart > 0);
        BOUT_FOR(i, g_22().getRegion("RGN_NOY")) {
          _g_22_yhigh.value()[i] =
              SQ(0.5 * (std::sqrt(g_22()[i]) + std::sqrt(g_22()[i.yp()])));
        }
      }
    }
  }
  return g_22_yhigh();
}

void Coordinates::_compute_cell_area_x() const {
  BOUT_OMP_SAFE(critical)
  {
    if (!_cell_area_xlow.has_value()) {
      const FieldMetric area_centre = J() / sqrt(g_11()) * dy_ * dz_;
      _cell_area_xlow.emplace(emptyFrom(area_centre));
      _cell_area_xhigh.emplace(emptyFrom(area_centre));
      // We cannot setLocation, as that would trigger the computation of staggered
      // metrics.
      ASSERT0(localmesh->xstart > 0);
      BOUT_FOR(i, area_centre.getRegion("RGN_NOX")) {
        (*_cell_area_xlow)[i] = 0.5 * (area_centre[i] + area_centre[i.xm()]);
        (*_cell_area_xhigh)[i] = 0.5 * (area_centre[i] + area_centre[i.xp()]);
      }
    }
  }
}

void Coordinates::_compute_cell_area_y() const {
  BOUT_OMP_SAFE(critical)
  {
    if (!_cell_area_ylow.has_value()) {
      if (g_11().isFci()) {
        const FieldMetric jxz_centre = J() / sqrt(g_22());
        auto jxz_ylow = emptyFrom(jxz_centre);
        auto jxz_yhigh = emptyFrom(jxz_centre);

        auto By_c = emptyFrom(jxz_centre);
        auto By_h = emptyFrom(jxz_yhigh);
        auto By_l = emptyFrom(jxz_ylow);
        if (localmesh->get(By_c, "By", 0.0, false, CELL_CENTRE) != 0) {
          throw BoutException("The grid file does not contain `By`.");
        }
        if (localmesh->get(By_l, "By_cell_ylow", 0.0, false) != 0) {
          throw BoutException("The grid file does not contain `By_cell_ylow`.");
        }
        if (localmesh->get(By_h, "By_cell_yhigh", 0.0, false) != 0) {
          throw BoutException("The grid file does not contain `By_cell_yhigh`.");
        }
        BOUT_FOR(i, By_c.getRegion("RGN_NOY")) {
          jxz_ylow[i] = By_c[i] / By_l[i] * jxz_centre[i];
          jxz_yhigh[i] = By_c[i] / By_h[i] * jxz_centre[i];
        }
        ASSERT3(isUniform(dx_, true, "RGN_ALL"));
        ASSERT2(isUniform(dx_, false, "RGN_ALL"));
        ASSERT3(isUniform(dz_, true, "RGN_ALL"));
        ASSERT2(isUniform(dz_, false, "RGN_ALL"));
        _cell_area_ylow.emplace(jxz_ylow * dx_ * dz_);
        _cell_area_yhigh.emplace(jxz_yhigh * dx_ * dz_);
      } else {
        // Field aligned
        const FieldMetric area_centre = J() / sqrt(g_22()) * dx_ * dz_;
        _cell_area_ylow.emplace(emptyFrom(area_centre));
        _cell_area_yhigh.emplace(emptyFrom(area_centre));
        // We cannot setLocation, as that would trigger the computation of staggered
        // metrics.
        BOUT_FOR(i, localmesh->getRegion("RGN_ALL")) {
          if (i.y() > 0) {
            (*_cell_area_ylow)[i] = 0.5 * (area_centre[i] + area_centre[i.ym()]);
          } else {
            (*_cell_area_ylow)[i] = BoutNaN;
          }
          if (i.y() < localmesh->LocalNy - 1) {
            (*_cell_area_yhigh)[i] = 0.5 * (area_centre[i] + area_centre[i.yp()]);
          } else {
            (*_cell_area_yhigh)[i] = BoutNaN;
          }
        }
      }
    }
  }
}

void Coordinates::_compute_cell_area_z() const {
  BOUT_OMP_SAFE(critical)
  {
    if (!_cell_area_zlow.has_value()) {
      const FieldMetric area_centre = J() / sqrt(g_33()) * dx_ * dy_;
      _cell_area_zlow.emplace(emptyFrom(area_centre));
      _cell_area_zhigh.emplace(emptyFrom(area_centre));
      // We cannot setLocation, as that would trigger the computation of staggered
      // metrics.
      BOUT_FOR(i, area_centre.getRegion("RGN_NOZ")) {
        (*_cell_area_zlow)[i] = 0.5 * (area_centre[i] + area_centre[i.zm()]);
        (*_cell_area_zhigh)[i] = 0.5 * (area_centre[i] + area_centre[i.zp()]);
      }
    }
  }
}

void Coordinates::_compute_cell_volume() const {
  BOUT_OMP_SAFE(critical)
  {
    if (!_cell_volume.has_value()) {
      _cell_volume.emplace(*jacobian_cache * dx_ * dy_ * dz_);
    }
  }
}

std::shared_ptr<YBoundary> Coordinates::makeYBoundary(YBndryType type) const {
  return std::make_shared<YBoundary>(type, localoptions, *localmesh);
}

void Coordinates::setBxy(FieldMetric Bxy, const bool communicate) {
  //TODO: Calculate Bxy and check value is close
  Bxy_ = std::move(Bxy);
  _g_22_ylow.reset();
  _g_22_yhigh.reset();
  invalidateAccessorCache();
  if (communicate) {
    localmesh->communicate_no_slices(Bxy_);
  }
}

void Coordinates::setContravariantMetricTensor(
    const ContravariantMetricTensor& metric_tensor, const std::string& region,
    bool recalculate_staggered, bool force_interpolate_from_centre) {
  contravariantMetricTensor = metric_tensor;
  covariantMetricTensor = contravariantMetricTensor.inverse(region);
  invalidateMetricCaches();
  setJ(recalculateJacobian());
  setBxy(recalculateBxy());
  recalculateAndReset(recalculate_staggered, force_interpolate_from_centre);
}

void Coordinates::setCovariantMetricTensor(const CovariantMetricTensor& metric_tensor,
                                           const std::string& region,
                                           bool recalculate_staggered,
                                           bool force_interpolate_from_centre) {
  covariantMetricTensor = metric_tensor;
  contravariantMetricTensor = covariantMetricTensor.inverse(region);
  invalidateMetricCaches();
  setJ(recalculateJacobian());
  setBxy(recalculateBxy());
  recalculateAndReset(recalculate_staggered, force_interpolate_from_centre);
}

void Coordinates::setMetricTensor(
    const ContravariantMetricTensor& contravariant_metric_tensor,
    const CovariantMetricTensor& covariant_metric_tensor) {
  contravariantMetricTensor = contravariant_metric_tensor;
  covariantMetricTensor = covariant_metric_tensor;
  invalidateMetricCaches();
  setJ(recalculateJacobian());
  setBxy(recalculateBxy());
}

void Coordinates::communicateMetricTensor() {
  contravariantMetricTensor.communicate();
  covariantMetricTensor.communicate();
}

void Coordinates::communicateDz() { localmesh->communicate(dz_); }

void Coordinates::splitBxyParallelSlices() {
  Bxy_.splitParallelSlices();
  Bxy_.yup() = Bxy_.ydown() = Bxy_;
}
